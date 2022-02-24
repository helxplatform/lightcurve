
import os, configparser, time, zlib
import boto3
import botocore
import numpy as np
import astropy_healpix
import astropy.units as u
import msgpack
import astropy.coordinates
import healpy
import os.path
import s3_util
import tempfile

local_config = configparser.ConfigParser()
local_config.read('config.ini')

class light_curve:
    def __init__(self, mjds, mags, magerrs, limmags, quality, cameras, ratchetids, files):
            self.mjds = mjds
            self.mags = mags
            self.magerrs = magerrs
            self.limmags = limmags
            self.quality = quality
            self.cameras = cameras
            self.ratchetids = ratchetids
            self.filenames = files

class s3_storage:
    def __init__(self):
        self.client = None
        self.last_client_use = 0
        self.ratchet_img_cache_loc = "/media/ram/ratchet_imgs/"
        self.s3aki, self.s3sak = self.get_k8s_s3secrets()

    def get_k8s_s3secrets(self):
        s3aki = os.environ.get('S3_ACCESS_KEY_ID', '')
        s3sak = os.environ.get('S3_SECRET_ACCESS_KEY', '')
        return (s3aki,s3sak)

    def get_client(self):
        if self.client is None or (time.time() - self.last_client_use) > 60.0:
            session = boto3.session.Session(
                aws_access_key_id=self.s3aki,
                aws_secret_access_key=self.s3sak
            )
            self.client = session.client(service_name='s3', endpoint_url='https://na-s3.renci.org')
        return self.client

    def get_s3resource(self):
        session = boto3.session.Session(
                aws_access_key_id=self.s3aki,
                aws_secret_access_key=self.s3sak
        )
        s3resource = session.resource(service_name='s3', endpoint_url='https://na-s3.renci.org')
        return s3resource

    def put_object(self, filename, obj_bytes):
        localfilename = "./tmp/s3upload.tmp"
        with open(localfilename, 'wb') as f:
            f.write(obj_bytes)
        try:
            s3client = self.get_client()
            s3client.upload_file(localfilename, "lightcurve", filename)
        except botocore.exceptions.ClientError as e:
            print(e.response)
            return
        return

    def delete_object(self, filename):
        s3client = self.get_client()
        try:
            s3client.delete_object(Bucket='lightcurve', Key=filename)
        except botocore.exceptions.ClientError as e:
            print(e.response)
            return
        return

    def download_object(self, filename):
        s3client = self.get_client()
        temp_name = next(tempfile._get_candidate_names()) # e.g. px9cp65s
        #localfilename = f"./tmp/s3obj-{temp_name}.tmp"
        localfilename = f"/tmp/s3obj-{temp_name}.tmp"
        try:
            # e.g. "light_curves/healpixes/000000.evrlc"
            # e.g. "light_curves/minipix2/000002.evrlc"
            # e.g. "light_curves/ratchet_imgs/ML0014214_20170128003309.jpg"
            s3filename = "light_curves/" + filename
            s3client.download_file('lightcurve', s3filename, localfilename)
        except botocore.exceptions.ClientError as e:
            if e.response['Error']['Code'] == "404":
                print(f"The object does not exist: {s3filename}")
            else:
                print("Some other error")
            return None
        with open(localfilename, 'rb') as f:
            s = f.read()
        os.remove(localfilename)
        return s

    def list_files(self):
        s3resource = self.get_s3resource()
        s3bucket = s3resource.Bucket('lightcurve')
        file_objects = s3bucket.objects.all()
        filenames = []
        #sizes = []
        #lastmodifieds = []
        for o in file_objects:
            filenames.append(o.key)
            #sizes.append(o.size)
            #lastmodifieds.append(o.last_modified)
        return filenames

    def upload_minipix(self, minipix, lc_bytes):
        self.put_object("minipix/%06i.evrlc" % (minipix), lc_bytes)

    def download_minipix(self, minipix):
        try:
            return self.download_object("minipix/%06i.evrlc" % (minipix))
        except google.api_core.exceptions.NotFound:
            return None

    def upload_minipix2(self, minipix, lc_bytes):
        self.put_object("minipix2/%06i.evrlc" % (minipix), lc_bytes)

    def download_minipix2(self, minipix):
        try:
            return self.download_object("minipix2/%06i.evrlc" % (minipix))
        except google.api_core.exceptions.NotFound:
            return None

    def download_healpix_metadata(self, healpix):
        try:
            return self.download_object("healpixes/%06i.evrlc" % (healpix))
        ###except google.api_core.exceptions.NotFound:
        except Exception as e:
            print(e)
            return None

    def download_epoch_map(self):
        return self.download_object("epoch_map.png")

    def download_ratchet_image(self, camera, ratchetid):
        return self.download_object("ratchet_imgs/ML%s_%s.jpg" % (camera, ratchetid))

    def download_ratchet_image_header(self, camera, ratchetid):
        return zlib.decompress(self.download_object("ratchet_imgs/ML%s_%s.header" % (camera, ratchetid)))

    def delete_files_containing_target(self, target):
        for fn in self.list_files():
            if target in fn:
                self.delete_object(fn)

    def unpack_light_curve(self,data,cameras,ratchetids,orig_epoch_n=None,filenames=None):
        minipix_data = np.frombuffer(data, dtype=np.uint16)
        data = np.reshape(minipix_data, (-1, 5))

        mjds = data[:,0] + ((data[:,1].astype(np.float64))/32767.0)
        mags = data[:,2].astype(np.float32) / 1000.0
        magerrs = (data[:,3] >> 8).astype(np.float32) / 255.0
        limmags = ((data[:,3] - ((data[:,3] >> 8) << 8)).astype(np.float32)/255.0)*10.0 + 10.0
        quality = data[:,4].astype(np.float32) / 255.0

        cameras = np.frombuffer(cameras,dtype=np.int64)
        ratchetids = np.frombuffer(ratchetids,dtype=np.int64)
        mask = (mags > 1.0) & (magerrs < 0.9)
        mjds = mjds[mask]
        mags = mags[mask]
        magerrs = np.sqrt((magerrs[mask] ** 2) + 0.01**2)
        limmags = limmags[mask]
        quality = quality[mask]

        orig_epoch_n = np.frombuffer(orig_epoch_n,dtype=np.int32)
        print (orig_epoch_n.shape, orig_epoch_n.dtype)
        files = filenames[orig_epoch_n]
        cameras = cameras[orig_epoch_n]
        ratchetids = ratchetids[orig_epoch_n]
        cameras = cameras[mask]
        ratchetids = ratchetids[mask]
        files = files[mask]

        if len(mjds) > 10:
            return light_curve(mjds, mags, magerrs, limmags, quality, cameras, ratchetids, files)
        else:
            return None

    def get_healpix_starlist(self, ra, dec):
        #, return_list_all=False, return_minipix_boundaries=False
        healpix = healpy.ang2pix(32, np.radians((dec*-1.0) + 90.0), np.radians(360.0 - ra))
        #healpix_conv = astropy_healpix.HEALPix(nside=32, order='ring')
        #healpix = healpix_conv.lonlat_to_healpix(ra * u.deg, dec * u.deg, return_offsets=False)


        healpix_metadata = self.download_healpix_metadata(healpix)

        if healpix_metadata is None:
            return None

        healpix_metadata = zlib.decompress(healpix_metadata)
        healpix_starprops, best_epoch, healpix_cameras, healpix_ratchetids, healpix_files = msgpack.unpackb(healpix_metadata,raw=False)
        healpix_files = np.frombuffer(healpix_files, dtype='<U62')
        healpix_cameras = np.frombuffer(healpix_cameras, dtype='int')
        healpix_ratchetids = np.frombuffer(healpix_ratchetids, dtype='int')

        star_list = []

        target = astropy.coordinates.SkyCoord(ra*u.deg, dec*u.deg, frame='fk5')

        for n,star in enumerate(healpix_starprops):
            coords = astropy.coordinates.SkyCoord(star[1]*u.deg, star[2]*u.deg, frame='fk5')
            sep = target.separation(coords)
            star_list.append([sep, n, star[1], star[2], star[3], star[7], star[5], star[6], star[4]])

        """
        if return_minipix_boundaries:
            bounds_ra, bounds_dec = healpix_conv.boundaries_lonlat([healpix], step=3)
            bounds_ra = bounds_ra.to(u.deg).value
            bounds_dec = bounds_dec.to(u.deg).value
            vertices = np.vstack([bounds_ra.ravel(), bounds_dec.ravel()]).transpose()
            return star_list, vertices
        """
        return star_list

    def get_light_curve(self, ra, dec, min_rad=1*u.arcmin):
        minipix_conv = astropy_healpix.HEALPix(nside=512, order='ring')
        minipix = minipix_conv.lonlat_to_healpix(ra * u.deg, dec * u.deg, return_offsets=False)
        healpix = healpy.ang2pix(32, np.radians((dec*-1.0) + 90.0), np.radians(360.0 - ra))
        #healpix_conv = astropy_healpix.HEALPix(nside=32, order='ring')
        #healpix = healpix_conv.lonlat_to_healpix(ra * u.deg, dec * u.deg, return_offsets=False)

        healpix_metadata = self.download_healpix_metadata(healpix)

        if healpix_metadata is None:
            return None

        healpix_metadata = zlib.decompress(healpix_metadata)

        healpix_starprops, best_epoch, healpix_cameras, healpix_ratchetids, healpix_files = msgpack.unpackb(healpix_metadata,raw=False)
        healpix_files = np.frombuffer(healpix_files, dtype='<U62')
        healpix_cameras = np.frombuffer(healpix_cameras, dtype='int')
        healpix_ratchetids = np.frombuffer(healpix_ratchetids, dtype='int')

        # should be all minipixes around the target point
        minipix_data = self.download_minipix2(minipix)
        if minipix_data is None:
                return None

        minipix_data = zlib.decompress(minipix_data)
        #       print (minipix_data, memoryview(minipix_data), memoryview(minipix_data).itemsize)
        stars = msgpack.unpackb(minipix_data,raw=False)[0]

        min_dist = 1e5*u.arcmin
        target_n = -1
        target_seps = []
        target = astropy.coordinates.SkyCoord(ra*u.deg, dec*u.deg, frame='fk5')

        for n,star in enumerate(stars):
            star_id = star[0]
            star_ra = star[1]
            star_dec = star[2]

            coords = astropy.coordinates.SkyCoord(star_ra*u.deg, star_dec*u.deg, frame='fk5')
            sep = target.separation(coords)
            target_seps.append(sep)

            if sep < min_dist:
                min_dist = sep
                target_n = n

        if min_dist < min_rad:
            return self.unpack_light_curve(data=stars[target_n][4],cameras=healpix_cameras,ratchetids=healpix_ratchetids,orig_epoch_n=stars[target_n][3],filenames=healpix_files)
        else:
            return None


if __name__ == "__main__":
    test = s3_storage()

    #test.list_files()
    #value = test.get_k8s_s3secrets()

    util = s3_util.s3_util()
    files, sizes, lastModifieds = util.list_folder_files('lightcurve', 'light_curves/ratchet_imgs', '.junk')
    test.put_object('light_curves/ratchet_imgs/blah.junk', bytes("stuff\n", 'UTF-8'))
    files, sizes, lastModifieds = util.list_folder_files('lightcurve', 'light_curves/ratchet_imgs', '.junk')
    test.delete_object('light_curves/ratchet_imgs/blah.junk')
    files, sizes, lastModifieds = util.list_folder_files('lightcurve', 'light_curves/ratchet_imgs', '.junk')

    """
    test.upload_star(3000, 1234, b"tgregweg")
    print (test.list_files())

    test.delete_all_stars_in_healpix(3000)
    print (test.list_files())

    test.put_object("test.txt", b"Hello world")
    """
    for l in open("ums_stars.tsv"):
        if len(l.strip()) > 0 and l[0] != '#':
            ra = float(l.split()[0])
            dec = float(l.split()[1])
            lc_out_fn = "lc_%03.4f_%02.4f.txt" % (ra,dec)
            if not os.path.exists(lc_out_fn):
                out_light_curve = test.get_light_curve(ra,dec)
                if out_light_curve is not None:
                    out = open("lc_%03.4f_%02.4f.txt" % (ra,dec), "w")
                    for mjd, mag, mag_err, quality in zip(out_light_curve.mjds, out_light_curve.mags, out_light_curve.magerrs, out_light_curve.quality):
                        print ("%.5f %.3f %.3f %f" %(mjd, mag, mag_err, quality), file=out)
                    print (ra,dec, len(out_light_curve.mjds))
                else:
                    print (ra, dec, "not found")
