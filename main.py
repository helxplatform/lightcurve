import datetime, json, io, pickle, os, time, sys, base64, csv

from flask import Flask, send_file

from PIL import Image, ImageDraw

import jinja2

from bokeh.embed import json_item
from bokeh.plotting import figure
from bokeh.layouts import column, row
from bokeh.resources import CDN
from bokeh.events import DoubleTap
from bokeh.models import CustomJS, ColumnDataSource, TapTool, Circle, RangeTool, WheelZoomTool, NumeralTickFormatter, Slider, TextInput, CheckboxButtonGroup
from bokeh.palettes import Spectral5
 
import astropy.wcs
import astropy.units as u
import astropy.coordinates

import s3_cloud
import logging
import numpy as np

app = Flask(__name__)

# test with python main.py
# Deploy to App Engine using: https://cloud.google.com/appengine/docs/standard/python3/building-app/deploying-web-service
# gcloud app deploy
# gcloud app browse
# log viewer page: https://console.cloud.google.com/logs/viewer
def targetid_from_radec(ra,dec):
    return "EVRJ%.3f%+.3f"%(ra,dec)

def jinga_render(template_name, **kwargs):
    templateLoader = jinja2.FileSystemLoader(searchpath="./templates/")
    templateEnv = jinja2.Environment(loader=templateLoader, autoescape=True) # autoescape off needed for header include, but may be vulnerable
    template = templateEnv.get_template(template_name)
    return template.render(kwargs)

def check_rate_limit():
    # this is a kill switch if something decides to crawl the entire DB
    # will only kill this one instance (although others should rapidly die if this is a real problem...)
    # /tmp is stored in a ramdisk per instance
    
    try:
        last_access_time = os.path.getmtime("/tmp/access_log.txt")
    except FileNotFoundError:
         last_access_time = time.time()

    if time.time() - last_access_time > 3600.0:
        with open("/tmp/access_log.txt","w") as last_access_log:
            print (datetime.datetime.now(), file=last_access_log)
    else:
        with open("/tmp/access_log.txt","a") as last_access_log:
            print (datetime.datetime.now(), file=last_access_log)

    with open("/tmp/access_log.txt","r") as last_access_log:
        n_recent_connections = len(last_access_log.readlines())
        if n_recent_connections > 40000:
            print ("Exceeded maximum connections, with %i connections" % n_recent_connections)
            sys.exit(1)

def sigma_clip(x,sigma,niter):
    x = np.array(x)
    if len(x) > 3:
       for i in range(niter):
           xt = x-np.mean(x)
           x = x[np.where(abs(xt) < sigma*np.std(xt))]
    return x



def make_lc_download(ra, dec, format):
    cloudstor = s3_cloud.s3_storage()

    lc = cloudstor.get_light_curve(ra, dec)

    mjds = lc.mjds
    mags = lc.mags 
    magerrs = lc.magerrs
    limmags = lc.limmags
    qualities = lc.quality 
    cameras = lc.cameras
    ratchets = lc.ratchetids

    targetid = targetid_from_radec(ra,dec)

    if format == "txt":
        out_string = io.StringIO()
        print ("#","%10s"%"MJD", "%8s"%"g-mag", "%8s"%"magerr", "%6s"%"limmag", "%5s"%"quality", "%7s"%"camera", "ratchetID", file=out_string)
        for epoch, mag, magerr, limmag, quality, camera, ratchet in zip(mjds,mags,magerrs,limmags,qualities,cameras,ratchets):
            print ("%12.6f"%epoch, "%8.4f"%mag, "%8.4f"%magerr, "%6.2f"%limmag, "%6.2f"%quality, "ML"+"%i"%camera, ratchet, file=out_string)

        # encode to bytes-format output
        mem = io.BytesIO()
        mem.write(out_string.getvalue().encode('utf-8'))
        mem.seek(0)
        return send_file(mem, as_attachment=True, download_name=targetid+".txt",mimetype='text/plain')

    if format == "csv":
        out_string = io.StringIO()
        writer = csv.writer(out_string)
        writer.writerow(["MJD", "g-mag", "magerr", "limmag", "quality", "camera", "ratchetID"])
        for epoch, mag, magerr, limmag, quality, camera, ratchet in zip(mjds,mags,magerrs,limmags,qualities,cameras,ratchets):
            writer.writerow([epoch, mag, magerr, limmag, quality, camera, ratchet])

        # encode to bytes-format output
        mem = io.BytesIO()
        mem.write(out_string.getvalue().encode('utf-8'))
        mem.seek(0)
        out_string.close()
        return send_file(mem, as_attachment=True, download_name=targetid+".csv",mimetype='text/csv')

    return None

def make_target_display_data(ra, dec, full_targets_table):
    print ("Making plot")
    t = time.time()

    cloudstor = s3_cloud.s3_storage()

    lc = cloudstor.get_light_curve(ra, dec)
    print ("Got LC", time.time() - t)

    mjds = lc.mjds
    mags = lc.mags 
    magerrs = lc.magerrs
    limmags = lc.limmags
    quality = lc.quality 
    cameras = lc.cameras
    ratchets = lc.ratchetids
    filenames = lc.filenames

    if mjds is None:
            mjds = np.zeros(1)
            mags = np.zeros(1)
            magerrs = np.zeros(1)
            limmags = np.zeros(1)
            quality = np.zeros(1)
            cameras = np.zeros(1)
            ratchets = np.zeros(1)
            filenames = np.zeros(1)
 
    mask = (magerrs < 1.0) & (mags > 5) & (quality < 10000.0)
    mjds = mjds[mask]
    mags = mags[mask]
    magerrs = magerrs[mask]
    quality = quality[mask]
    limmags = limmags[mask]
    cameras = cameras[mask]
    ratchets = ratchets[mask]
    filenames = filenames[mask]

    #print ("Mageers:", magerrs, "libmags:", limmags)

    best_ratchetcamera_n = np.argmin(magerrs)
    best_ratchetcamera = "%i_%i" % (cameras[best_ratchetcamera_n], ratchets[best_ratchetcamera_n])

    cameraratchets = []
    for c,r in zip(cameras,ratchets):
        cameraratchets.append("%i_%i"%(c,r))

    mags_middle = np.median(mags)
    mags_high = np.max(mags) - mags_middle
    mags_low = mags_middle - np.min(mags)
    if mags_high < 0.5:
        mags_high = 0.5
    if mags_low < 0.5:
        mags_low = 0.5

    p = figure(title = "", 
                sizing_mode="scale_both", 
                plot_width=1000,
                plot_height=300,
                toolbar_location="right",
                tools="box_zoom, pan, reset",
                active_drag="box_zoom",
                #x_axis_location="above",
                background_fill_color="#efefef",
                x_range=(np.min(mjds)-10,np.max(mjds)+10),
                y_range=(mags_middle + mags_high*1.5, mags_middle - mags_low*1.5),
                margin=(10,10,10,10),
                output_backend="webgl"
                )

    p.axis.axis_label_text_font_style = 'bold'
    p.axis.axis_label_text_font_size = "12pt"
    p.axis.axis_label_text_font = "\"Helvetica Neue\", \"Open Sans\", Helvetica, Arial, sans-serif"
    
    p.background_fill_color = "#F8F8FF"
    p.yaxis[0].formatter = NumeralTickFormatter(format="0.00")

    p.border_fill_alpha = 0.0
    p.outline_line_color = "black"

    p.toolbar.active_scroll = p.select_one(WheelZoomTool)

    p.xaxis.axis_label = "MJD"
    p.yaxis.axis_label = "Evryscope-mag"
    p.yaxis.major_label_text_font_size = "12pt"
    p.xaxis.major_label_text_font_size = "12pt"
    p.toolbar.logo = None


    good_color_r = 0.0
    good_color_g = 162.0/255.0
    good_color_b = 252.0/255.0

    color_delta = (quality - 0.4) / 0.6

    color_delta[quality < 0.4] = 0.0

    color_r = good_color_r * (1.0 - color_delta) + 1.0 * color_delta
    color_g = good_color_g * (1.0 - color_delta) + 0.9 * color_delta
    color_b = good_color_b * (1.0 - color_delta) + 0.85 * color_delta
    
    data = {
    'x': mjds,
    'y': mags,
    'mjds': mjds,
    'mags': mags,
    'quality': quality,
    'cameraratchets': cameraratchets,
    'color': ["#%02x%02x%02x" % (int(r_color*255.0), int(g_color*255.0), int(b_color*255.0)) for r_color, g_color, b_color in zip(color_r, color_g, color_b)]
    #'star_ra': np.ones(len(mjds))*ra,
    #'star_dec': np.ones(len(mjds))*dec # for later position changes with time
    }
    
    target_info = ""
    target_info+= "<table class=\"targettable\">"
    target_info+= "<tr><td><b>Position:</b></td> <td>%.4f, %.4f</td></tr>" % (ra,dec)
    target_info+= "<tr><td><b>Magnitude:</b></td> <td>%.2f (Evryscope-g)</td></tr>" % (np.median(mags))
    target_info+= f"<tr><td><b>Epochs:</b></td> <td>{len(mjds):,} 2-min exposures</td></tr>"
    target_info+= "<tr><td><b>Survey length: </b></td> <td>%i d</td></tr>" % int((np.max(mjds) - np.min(mjds)))
    target_info+= "<tr><td><b>Nights with data:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </b></td> <td>%i d</td></tr>" % len(set(mjds.astype(np.int32)))
    target_info+= "</table>"
    data = ColumnDataSource(data)

    click_point_callback = CustomJS(args={"current_src": data}, code=open("templates/get_image.js","r").read())
    tap = TapTool(callback=click_point_callback)

    p.circle('x', 'y', size=5.0, fill_alpha=0.5, source=data, line_color=None,color='color',nonselection_fill_color='color')
    p.tools.append(tap)

    #period_slider = Slider(start = -0.01, end = 0.01, value=0.0, step=0.00001, title=None, disabled=True, width=120, format="0.00000f")
    bin_slider = Slider(start = 5, end = 100, value=5, step=5, title="Epoch binning", disabled=True, width=120, format="i")

    period_text = TextInput(placeholder='Phasing period (d)', disabled=True, width=150)

    button_group = CheckboxButtonGroup(labels=["Show only higher-quality points", "Phase light curve", "Bin points"],width=500)
    
    lc_plot_update_js = """
    Math.fmod = function (a,b) { return Number((a - (Math.floor(a / b) * b))) };
    var period_fold_on = false;
    var bin_on = false;
    var only_highq_points = false;

    var selections = Array.from(buttons.active);
    //console.log(typeof(buttons.active), typeof(buttons.active[0]), Array.from(buttons.active));
    
    if (selections.includes(0)) only_highq_points = true;
    if (selections.includes(1)) period_fold_on = true;
    if (selections.includes(2)) bin_on = true;

    //console.log(only_highq_points, bin_on, period_fold_on, this.active);

    if (period_fold_on) 
    {
       period_text.disabled = false;
       //period_slider.disabled = false; 
    }
    else
    {
       period_text.disabled = true;
       //period_slider.disabled = true;
    }

    if (bin_on == true) bin_slider.disabled = false;
    else bin_slider.disabled = true;

    // make the points for the graph depending on the settings
    // (update the graph regardless of if necessary; inefficient)
    const x = lc_data.data['x'];
    const y = lc_data.data['y'];
    const epochs = lc_data.data['mjds'];
    const quality = lc_data.data['quality'];
    const mags = lc_data.data['mags'];
    const period = parseFloat(period_text.value);

    if (period_fold_on == false)
    {
           if (plotRange.start < 10000.0)
           {
              plotRange.start = epochs[0];
              plotRange.end = epochs[epochs.length - 1];
           }

           for (var i=0; i < x.length; i++) 
           {
               x[i] = epochs[i];
               y[i] = mags[i];
           }
    }
    else
    {
         if (period > 0.0)
         {
             for (var i=0; i < x.length; i++)
             {
                x[i] = Math.fmod(epochs[i], period);
                y[i] = mags[i];
             }
             if (plotRange.start > 40000.0)
             {
                 plotRange.start = 0.0;
                 plotRange.end = period;
             }
             //const period_slider_mid = Math.floor(period*1000.0)
             //period_slider.start = (period_slider_mid - 0.49)/1000.0;
             //period_slider.end = (period_slider_mid + 0.49)/1000.0;
             //period_slider.value = period;
             //period_slider.step = period*1e-8;
         }           
    }

    if (only_highq_points == true)
    { 
        for (var i=0; i < y.length; i++) 
        {
            if (quality[i] > 0.4) y[i] = 1/0;
        }
    }

    if (bin_on == true)
    {
        if (period_fold_on == true)
        {
            var array_x = Array.from(x);
            var array_y = Array.from(y);
            let argsort = a=>a.map(d).sort().map(u);let d=(v,i)=>[v,i];let u=i=>i[1];
            let order = argsort(array_x);            
            array_x = order.map(i => x[i]);
            array_y = order.map(i => y[i]);
            for (var i=0; i < array_x.length; i++)
            {
                x[i] = array_x[i];
                y[i] = array_y[i];
            }
        }

        const bin_n = bin_slider.value;
        var avg = 0.0;
        var n_points = 0;
        var n_in_bin = 0;
        for (var i=0; i < x.length; i++)
        {
           if (y[i] > 5.0 && y[i] < 25.0) 
           {
              n_in_bin+=1;
              avg+=y[i];
           }
           n_points+=1;
           
           if (n_points >= bin_n)
           {
               y[Math.floor(i / bin_n)*bin_n] = avg / n_in_bin;
               for (var i2=1;i2<bin_n;i2++) y[Math.floor(i/bin_n)*bin_n + i2] = 1/0;
               n_points = 0;
               avg = 0.0;  
               n_in_bin = 0;
           }
        }

        for (var i2=Math.floor(i/bin_n)*bin_n;i2<y.length;i2++) y[i2] = 1/0;
    }

    lc_data.change.emit();
    """

    button_group.js_on_click(CustomJS(args=dict(lc_data=data, period_text=period_text, bin_slider = bin_slider, plotRange=p.x_range, buttons=button_group), code=lc_plot_update_js))

    #period_slider.js_on_change('value', CustomJS(args=dict(data=data, period_text=period_text, period_slider=period_slider,plotRange=p.x_range), code=lc_plot_update_js))
    period_text.js_on_change('value', CustomJS(args=dict(lc_data=data, period_text=period_text, bin_slider=bin_slider, plotRange=p.x_range, buttons=button_group), code=lc_plot_update_js))
    bin_slider.js_on_change('value', CustomJS(args=dict(lc_data=data, period_text=period_text, bin_slider=bin_slider, plotRange=p.x_range, buttons=button_group), code=lc_plot_update_js))

    p.js_on_event(DoubleTap, CustomJS(args=dict(lc_data=data, period_text=period_text, plotRange_x=p.x_range, plotRange_y = p.y_range, buttons = button_group, p=p), code="""
    var period_fold_on = false;
    var bin_on = false;
    var only_highq_points = false;

    var selections = Array.from(buttons.active);
    
    if (selections.includes(0)) only_highq_points = true;
    if (selections.includes(1)) period_fold_on = true;
    if (selections.includes(2)) bin_on = true;
 
    const epochs = lc_data.data['mjds']; 
    const mags = lc_data.data['mags']; 
    const period = parseFloat(period_text.value);
   
    if (period_fold_on == true) {plotRange_x.start = 0.0; plotRange_x.end = period;}
    else {plotRange_x.start = epochs[0]; plotRange_x.end = epochs[epochs.length - 1];}    

    const filtered_mags = mags.filter( value => !Number.isNaN(value) );
    plotRange_y.end = Math.min(...filtered_mags)-0.5;
    plotRange_y.start = Math.max(...filtered_mags)+0.5; 
    console.log(plotRange_x.start, plotRange_x.end, plotRange_y.start, plotRange_y.end);
    """))    

    print ("Making table", time.time() - t)
    # make the light curve table
    lc_table = """
        <div class="scrollingtable" id="lc_table_scroll">
        <div id="lc_table_scroll2">
            <div id="lc_table_scroll3">
                <table id="lc_table_id">
                    <!-- <caption>Top Caption</caption> -->
                    <thead>
                        <tr>
                         <th><div label="MJD"></div></th>
                         <th><div label="g-mag"></div></th>
                         <th><div label="Magerr"></div></th>
                         <th><div label="Lim-mag"></div></th>
                         <th><div label="Quality (0..1)"></div></th>
                         <th><div label="Camera"></div></th>
                         <th><div label="RatchetID"></div></th>
                         <th><div label="Filename"></div></th>
              <th class="scrollbarhead"></th> <!--ALWAYS ADD THIS EXTRA CELL AT END OF HEADER ROW-->
                        </tr>
                    </thead>
                    <tbody>                        
    """

    n_dps = len(mjds)
    if full_targets_table == False:
        n_dps = 250

    for n in range(n_dps):
        lc_table+= "<tr><td>%.6f</td> <td>%.3f</td> <td>%.3f</td> <td>%.2f</td> <td>%.2f</td> <td>ML%s</td> <td>%s</td> <td>%s</td></tr>" % (mjds[n], mags[n], magerrs[n], limmags[n], quality[n], cameras[n], ratchets[n], filenames[n])
    
    lc_table+= """
               </tbody>
                </table>
    """
    
    if full_targets_table == False:
        lc_table+="<a href=\"javascript:void(0);\" onclick=\"update_table_to_full();\">Load the full dataset into table (another %i points) </a>" % (len(mjds) - n_dps)

    lc_table+="""
            </div>
        </div>
        </div>
    """

    print ("Finished", time.time() - t)

    return column(p,row(button_group,period_text,bin_slider)), best_ratchetcamera, target_info, lc_table

def make_image(ra, dec, cameraratchet, width, height, plot_stars):
    cloudstor = s3_cloud.s3_storage()
    camera, ratchetid = cameraratchet.split("_")
    
    image_data = cloudstor.download_ratchet_image(camera, ratchetid)
    header = cloudstor.download_ratchet_image_header(camera, ratchetid)
    wcs = astropy.wcs.WCS(header)

    tx, ty = wcs.all_world2pix(ra,dec,0)

    # make a large cutout around that area, for rotation
    image = Image.open(io.BytesIO(image_data))
    width = int(width)
    height = int(height)
    image = image.crop((tx-width/2,ty-height/2,tx+width/2,ty+height/2))
    image = image.convert("RGB")

    draw = ImageDraw.Draw(image)
    cx, cy = image.width / 2, image.height / 2
    marker_len = 10
    marker_skip = 5
    
    if plot_stars == True:
        # plot all nearby stars in DB
        cloudstor = s3_cloud.s3_storage()
        lcs = cloudstor.get_healpix_starlist(ra, dec)

        #mpix_boundaries_pixels = wcs.all_world2pix(mpix_boundaries,0)
        #mpix_boundaries_pixels[:,0]-=(tx-width/2)
        #mpix_boundaries_pixels[:,1]-=(ty-height/2)

        #for n in range(len(mpix_boundaries_pixels)-1):
        #    draw.line((mpix_boundaries_pixels[n,0],mpix_boundaries_pixels[n,1],mpix_boundaries_pixels[n+1,0],mpix_boundaries_pixels[n+1,1]), fill=(255,255,255),width=2)

        #draw.line((mpix_boundaries_pixels[0,0],mpix_boundaries_pixels[0,1],mpix_boundaries_pixels[len(mpix_boundaries_pixels)-1,0],mpix_boundaries_pixels[len(mpix_boundaries_pixels)-1,1]), fill=(255,255,255),width=2)

        for lc in lcs:
            lx, ly = wcs.all_world2pix(lc[2],lc[3],0)
            lx-=tx
            ly-=ty
            lx+=width/2
            ly+=height/2
            marker_rad = 4
            
            draw.ellipse((lx-marker_rad,ly-marker_rad,lx+marker_rad,ly+marker_rad), outline=(50,255,70))

    draw.line((cx - marker_skip, cy, cx - marker_skip - marker_len, cy), fill=(255,50,100),width=2)
    draw.line((cx + marker_skip, cy, cx + marker_skip + marker_len, cy), fill=(255,50,100),width=2)
    draw.line((cx, cy - marker_skip, cx, cy - marker_skip - marker_len), fill=(255,50,100),width=2)
    draw.line((cx, cy + marker_skip, cx, cy + marker_skip + marker_len), fill=(255,50,100),width=2)


    image_byte_array = io.BytesIO()
    image.save(image_byte_array,format='png',quality=99)

    return image_byte_array.getvalue()

def make_targets_table(ra,dec):
    cloudstor = s3_cloud.s3_storage()

    lcs = cloudstor.get_healpix_starlist(ra, dec)

    if lcs is not None and len(lcs) > 0:
        sidebar_items = ""

        out = """<div class="table_container"><div class="table">"""
        out+= """<div class="table-content">\n"""

        out+="""
                <div class="table-row-header">
                    <div class="table-data header-item"><a id="sep" class="filter__link filter__link--number" href="#">Sep. / arcsec</a></div>
                    <div class="table-data header-item"><a id="target_id" class="filter__link filter__link--number" href="#">ID</a></div>
                    <div class="table-data header-item"><a id="RA" class="filter__link filter__link--number" href="#">RA / deg</a></div>
                    <div class="table-data header-item"><a id="Dec" class="filter__link filter__link--number" href="#">Dec / deg</a></div>
                    <div class="table-data header-item"><a id="Evrymag" class="filter__link filter__link--number" href="#">Evry-mag (g)</a></div>
                    <div class="table-data header-item"><a id="Epochs" class="filter__link filter__link--number" href="#"># Epochs</a></div>
                    <div class="table-data header-item"><a id="StartMJD" class="filter__link filter__link--number" href="#">StartMJD</a></div>
                    <div class="table-data header-item"><a id="EndMJD" class="filter__link filter__link--number" href="#">EndMJD</a></div>
                    <div class="table-data header-item"><a id="magstd" class="filter__link filter__link--number" href="#">Mags-std</a></div>
                </div>
        """

        for lc in sorted(lcs):
            out+="""<div class="table-row">\n"""
            out+="""<div class="table-data">""" + ("%.2f" % lc[0].to_value(u.arcsec)) + "</div>\n"
            target_id = targetid_from_radec(lc[2],lc[3])
            out+="""<div class="table-data">""" + "<a href=\"\\dynamic\\target_display\\" + ("%.5f"%lc[2]) + "\\" + ("%.5f"%lc[3]) + ("\">%s</a>"%target_id)  + "</div>\n"
            out+="""<div class="table-data">""" + ("%.4f"%lc[2]) + "</div>\n"
            out+="""<div class="table-data">""" + ("%.4f"%lc[3]) + "</div>\n"
            out+="""<div class="table-data">""" + ("%.2f"%lc[4]) + "</div>\n"
            out+="""<div class="table-data">""" + ("%i"%lc[5]) + "</div>\n"
            out+="""<div class="table-data">""" + ("%.2f"%lc[6]) + "</div>\n"
            out+="""<div class="table-data">""" + ("%.2f"%lc[7]) + "</div>\n"
            out+="""<div class="table-data">""" + ("%.2f"%lc[8]) + "</div>\n"
            out+="</div>\n"
        out+="</div>\n"
        out+="</div></div>\n"



        sidebar_items+="<h3 class=\"sectionheading\"><span>Nearby stars list</span></h2>"
        sidebar_items+="<table class=\"targettable\"><tr><td>Query target:</td><td> %.4f, %.4f</td></tr>" % (ra, dec)
        galactic_coords = astropy.coordinates.SkyCoord(ra*u.deg, dec*u.deg, frame='icrs').galactic
        sidebar_items+="<tr><td>Galactic coordinates:</td><td> %.4f, %.4f</td></tr>" % (galactic_coords.l.degree, galactic_coords.b.degree)

        print (ra,dec)
        best_match_lc = cloudstor.get_light_curve(ra, dec, min_rad=1e5 * u.arcmin)

        # find best cameraratchet for closest match
        best_ratchetcamera_n = np.argmin(best_match_lc.magerrs)
        best_ratchetcamera = "%i_%i" % (best_match_lc.cameras[best_ratchetcamera_n], best_match_lc.ratchetids[best_ratchetcamera_n])

        sidebar_items+="<tr><td>Stars in closest minipix:</td><td> %i</td></tr>" % (len(lcs))
        sidebar_items+="</table>"
        sidebar_items+="<h3 class=\"sectionheading\" style=\"margin-top:25px\"><span>Evryscope image of region</span></h3>"

        sidebar_items+="<img width=\"100%\" src=\"/dynamic/image/" + repr(ra) + "/" + repr(dec) + "/" + best_ratchetcamera + "/300/200/plot_stars=True\">"
        sidebar_items+="<p class=\"expln_text\">red crosshair = ra,dec</p>"
        sidebar_items+="<p class=\"expln_text\">green circles = stars in DB</p>"
        sidebar_items+="<p class=\"expln_text\">white lines = boundary of query region</p>"

    else:
        out = "<div id=\"no_target_err\">No targets found</div>"
        sidebar_items = ""

    return out, sidebar_items

@app.route('/',methods=['GET'])
def root():
    check_rate_limit()
    return jinga_render("main_page.html", title="Evryscope Light Curves Server")

@app.route('/dynamic/target_list/<string:ra>/<string:dec>',methods=['GET'])
def target_list(ra,dec):
    print(f"target_list: {ra} {dec}", flush=True)
    check_rate_limit()
    ra = float(ra)
    dec = float(dec)
    target_list, sidebar_items = make_targets_table(ra,dec)
    return jinga_render("target_list.html", title="Evryscope targets for %.4f, %.4f" % (ra,dec), sidebar_items = sidebar_items, target_list = target_list)

@app.route('/dynamic/target_display/<string:ra>/<string:dec>',methods=['GET'])
def target_display(ra=None,dec=None):
    print(f"target_display: {ra} {dec}", flush=True)
    check_rate_limit()
    ra = float(ra)
    dec = float(dec)
    return jinga_render("target_display.html",resources=CDN.render(),star_ra=ra,star_dec=dec,sidebar_items="")
 
@app.route('/dynamic/target_display_data/<string:ra>/<string:dec>/<string:full>')
def target_display_data(ra,dec,full):
    print(f"target_display_data: {ra} {dec} {full}", flush=True)
    get_full_table_data = True
    if full != "full":
        get_full_table_data = False

    check_rate_limit()
    ra = float(ra)
    dec = float(dec)
    p, image_spec, target_info, lc_table = make_target_display_data(ra,dec,get_full_table_data)
    data = json.dumps([json_item(p, "lc_plot"), [ra,dec,image_spec], target_info, lc_table])
    print ("Data length:", len(data))
    return data

@app.route('/dynamic/download_lc/<string:ra>/<string:dec>/<string:format>')
def download_lc(ra,dec,format):
    check_rate_limit()
    ra = float(ra)
    dec = float(dec)
    output_format = None
    if format == "txt":
        output_format = "txt"
    elif format == "csv":
        output_format = "csv"
    elif format == "fits":
        output_format = "fits"

    if output_format is not None:
            return make_lc_download(ra,dec,format)
    else:
        return False

@app.route('/dynamic/image/<float(signed=True):ra>/<float(signed=True):dec>/<string:cameraratchet>/<int:width>/<int:height>/<string:plot_stars>/')
def image(ra,dec,cameraratchet,width,height,plot_stars): 
    print(f"image: {ra} {dec} {cameraratchet} {width} {height} {plot_stars}", flush=True)
    check_rate_limit()
    if plot_stars == "plot_stars=True":
        plot_stars_flag = True
    else:
        plot_stars_flag = False

    image = make_image(ra,dec,cameraratchet,width,height,plot_stars_flag)

    return send_file(io.BytesIO(image),
                     download_name='image.jpg',
                     mimetype='image/jpg')

@app.route('/dynamic/epoch_map')
def epoch_map(): 
    check_rate_limit()
    image = s3_cloud.s3_storage().download_epoch_map()
    
    return send_file(io.BytesIO(image),
                     download_name='image.png',
                     mimetype='image/png')

@app.route('/favicon.ico')
def favicon_ico(): 
    check_rate_limit()
    image = open("favicon.ico","rb").read()
    
    return send_file(io.BytesIO(image),
                     download_name='favicon.ico',
                     mimetype='image/ico')


if __name__ == '__main__':
    #cloudenv = os.environ.get('cloudenv', '')
    #print(f"cloudenv: '{cloudenv}'", flush=True)
    #if cloudenv == 'renci':

    print(f"app.static_url_path: '{app.static_url_path}'", flush=True)
    app.logger.setLevel(logging.DEBUG)
    app.run(host='0.0.0.0', port=8080, debug=True)

    #else:
    #    app.run(host='127.0.0.1', port=8080, debug=True)
