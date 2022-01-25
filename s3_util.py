
import boto3

class s3_util:
    def __init__(self):
        self.client = None

    def get_client(self):
        if self.client is None:
            session = boto3.session.Session(
                aws_access_key_id='',
                aws_secret_access_key=''
            )
            self.client = session.client(service_name='s3', endpoint_url='https://na-s3.renci.org')
        return self.client

    def list_folder_files(self, bucket:str, folder:str, extension:str):
        files = []
        sizes = []
        lastModifieds = []

        s3client = self.get_client()
        paginator = s3client.get_paginator('list_objects_v2')
        try:
            pages = paginator.paginate(Bucket=bucket, Prefix=folder)
            for page in pages:
                for obj in page['Contents']:
                    filename = obj['Key']
                    keep = True
                    if len(extension) > 0 and extension != '.*' and extension != '*':
                        if not filename.endswith(extension):
                            keep = False
                    if keep:
                        size = obj['Size']
                        lastModified = obj['LastModified']
                        files.append(filename)
                        sizes.append(size)
                        lastModifieds.append(lastModified)
        except botocore.exceptions.ClientError as e:
            print(f"response: {e.response}")
            return None
        return files, sizes, lastModifieds


if __name__ == "__main__":
    test = s3_util()
    files, sizes, lastModifieds = test.list_folder_files('lightcurve', 'light_curves/minipix', '.evrlc')
    print(f"len(files): {len(files)}, len(sizes): {len(sizes)}, len(lastModifides): {len(lastModifieds)}")
