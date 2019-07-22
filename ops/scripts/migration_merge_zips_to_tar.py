from typing import Dict
import zipfile
import tarfile
import os.path
import sys

from nomad import config


names: Dict[str, bool] = dict()


def add_upload(tf, upload):
    files = [
        os.path.join(config.fs.public, upload[0:2], upload, base)
        for base in ['raw-public.plain.zip', 'raw-restricted.plain.zip']]

    for f in files:
        with zipfile.ZipFile(f) as zf:
            for zipinfo in zf.infolist():
                name = zipinfo.filename
                if name not in names:
                    names[name] = True
                    with zf.open(zipinfo) as bf:
                        tarinfo = tarfile.TarInfo(name)
                        tarinfo.size = zipinfo.file_size
                        tf.addfile(tarinfo, bf)


target = sys.argv[1]

with tarfile.TarFile(target, 'x') as tf:
    for upload in sys.argv[2:]:
        print('adding upload %s' % upload)
        add_upload(tf, upload)
