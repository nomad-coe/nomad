#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

'''
Automatically synchronizes nomad it with a given database. It creates a list of paths
to mainfiles in nomad and compares it with paths in the external database. The missing
paths in nomad will then be downloaded from the external database and subsequently
uploaded to nomad. The downloaded files are by default saved in '/nomad/fairdi/external'.
'''

import requests
import re
import subprocess
from urllib import parse as urllib_parse
import os
import datetime
import click
import tarfile
import threading
import time
import typing

from .client import client
from nomad import config
from nomad.cli.client import upload as nomad_upload


class DbUpdater:
    def __init__(self, *args, **kwargs):
        self.db_name = 'aflowlib'
        self.root_url = 'http://aflowlib.duke.edu/AFLOWDATA/LIB1_LIB'
        self.local_path = '/nomad/fairdi/external'
        self.dbfile = None
        self.nomadfile = None
        self.outfile = None
        self.do_download = False
        self.do_upload = False
        self.do_publish = False
        self.cleanup = False
        self.parallel = 2
        self.max_depth = 10
        self.uids = []
        self._set(**kwargs)
        self._set_db()
        self._set_local_path()
        self.configure_client()

    def _set(self, **kwargs):
        for key, val in kwargs.items():
            key = key.lower()
            if hasattr(self, key):
                setattr(self, key, val)
            else:
                raise KeyError('Invalid key %s' % key)

    def configure_client(self, **kwargs):
        from nomad.cli.client import create_client
        self.client = create_client()

    def _set_local_path(self):
        subdir = ''
        if self.db_name.lower() == 'aflowlib':
            subdir = self.root_url.strip('/').split('/')[-1]
        dbpath = os.path.join(self.local_path, self.db_name)
        self._local_path = os.path.join(dbpath, subdir)
        if not os.path.isdir(dbpath):
            os.mkdir(dbpath)
        if not os.path.isdir(self._local_path):
            os.mkdir(self._local_path)

    def _set_db(self):
        if self.db_name.lower() == 'aflowlib':
            self.max_depth = 4
        else:
            raise NotImplementedError('%s not yet supported.' % self.db_name)

    def _open_page(self, path: str):
        with requests.Session() as session:
            response = session.get(path, verify=False)
        if not response.ok:
            return response.raise_for_status()
        return response

    def get_paths(self, root: str) -> typing.List[str]:
        response = self._open_page(root)
        paths = []
        for url in re.findall('<a href="([^"]+)">', str(response.content)):
            if url in response.url:
                continue
            paths.append(os.path.join(response.url, url.lstrip('./')))
        return paths

    def _is_mainfile(self, path: str) -> bool:
        if 'vasprun.xml' in path:
            return True
        return False

    def _is_dir(self, path: str) -> bool:
        path = path.strip()
        if path[-1] == '/' and self.root_url in path:
            return True
        return False

    def _depth(self, path: str) -> int:
        return urllib_parse.urlparse(path).path.rstrip('/').count('/')

    def _rules_ok(self, path: str) -> bool:
        ok = True
        if '?' in path:
            ok = False
        return ok

    def _to_namesafe(self, path: str) -> str:
        path = path.strip('/')
        return re.sub(r'[^\w\d-]', '_', path)

    def _read_from_file(self, filename):
        print('Reading from file %s' % filename)
        data = []
        with open(filename) as f:
            line = f.readline()
            while line:
                data_i = line.strip().split()

                if len(data_i) == 1:
                    data.append(data_i[0])

                else:
                    info = data_i[1] in ['T', 'True', '1', 'Y']
                    data.append((data_i[0], info))

                line = f.readline()
        return data

    def _write_to_file(self, data: typing.List, filename: str):
        with open(filename, 'w') as f:
            for i in range(len(data)):
                if isinstance(data[i], str):
                    f.write('%s\n' % data[i])
                else:
                    f.write('%s %s \n' % (data[i][0], data[i][1]))

    def _gen_db_list(self):
        print('Generating list from %s' % self.root_url)
        self.db_files = []
        todo = self.get_paths(self.root_url)
        while len(todo) > 0:
            cur = todo.pop(0)
            if not self._is_dir(cur):
                continue

            if self._depth(cur) == self.max_depth:
                self.db_files.append(cur)

            else:
                add = self.get_paths(cur)
                add = [path for path in add if self._rules_ok(path)]
                todo = add + todo

        if self.dbfile is not None:
            self._write_to_file(self.db_files, self.dbfile)

    def get_db_list(self):
        if self.dbfile is not None and os.path.isfile(self.dbfile):
            self.db_files = self._read_from_file(self.dbfile)
        else:
            self._gen_db_list()

    def _gen_nomad_list(self):
        print('Generating NOMAD list')
        if self.db_name.lower() == 'aflowlib':
            servers = ['LIB%d_LIB' % n for n in range(1, 10)]
            paths = []
            for s in servers:
                if s in self.root_url:
                    paths.append(s)
            if len(paths) == 0:
                paths = servers
            # uploader: Stefano Curtarolo
            kwargs = dict(
                uploader_id='81b96683-7170-49d7-8c4e-e9f34906b3ea',
                paths=paths,
                scroll=True)

        self.nomad_files = []
        while True:
            res = self.client.repo.search(**kwargs).response()
            results = res.result.results
            if len(results) == 0:
                break
            for i in range(len(results)):
                self.nomad_files.append(results[i]['mainfile'])
            scroll_id = res.result.scroll.scroll_id
            kwargs.update({'scroll_id': scroll_id})

        if self.nomadfile is not None:
            self._write_to_file(self.nomad_files, self.nomadfile)

    def get_nomad_list(self):
        if self.nomadfile is not None and os.path.isfile(self.nomadfile):
            self.nomad_files = self._read_from_file(self.nomadfile)
        else:
            self._gen_nomad_list()

    def _reduce_list(self, ilist: typing.List[str]):
        olist = []
        for e in ilist:
            p = urllib_parse.urlparse(e).path.strip('/')
            olist.append(os.path.join(*p.split('/')[1:self.max_depth]))
        olist = list(set(olist))
        olist.sort()
        return olist

    def compare_lists(self):
        print('Identifying differences')
        db = self._reduce_list(self.db_files)
        nomad = self._reduce_list(self.nomad_files)
        ns = set(nomad)
        self.update_list = [i for i in db if i not in ns]

        ds = set(db)
        self.in_nomad = [i for i in nomad if i not in ds]
        if len(self.in_nomad) > 0:
            fn = 'in_nomad.txt'
            print('Warning: Some NOMAD entries not found in db.')
            print('See %s for list.' % fn)
            self._write_to_file(self.in_nomad, fn)

        # add the root back
        u = urllib_parse.urlparse(self.root_url)
        up = u.path.strip('/').split('/')[0]
        root = '%s://%s/%s' % (u.scheme, u.netloc, up)
        self.update_list = [os.path.join(root, e) for e in self.update_list]
        self.is_updated_list = [False] * len(self.update_list)
        print('Found %d entries to be added in NOMAD' % len(self.update_list))

        if self.outfile is not None:
            data = [self.update_list[i] for i in range(len(self.update_list))]
            self._write_to_file(data, self.outfile)

    def _download(self, path: str, iodir: str) -> float:
        files = self.get_paths(path)
        files = [f for f in files if self._rules_ok(f)]
        size = 0.0
        with requests.Session() as session:
            for f in files:
                if self._is_dir(f):
                    _, s = self.get_files(f)
                    size += s
                    continue

                res = session.get(f, stream=True)
                fn = res.url.split('/')[-1]
                fn = os.path.join(iodir, fn)
                with open(fn, 'wb') as fb:
                    fb.write(res.content)

                size += os.path.getsize(fn)
        return size

    def _download_size(self, dirname: str) -> float:
        complete = False
        files = os.listdir(dirname)

        if self.db_name.lower() == 'aflowlib':
            if 'vasprun.xml.relax2.xz' in files:
                complete = True
        else:
            complete = True

        size = 0.0
        if not complete:
            self._cleanup([os.path.join(dirname, f) for f in files])
        else:
            for f in files:
                size += os.path.getsize(os.path.join(dirname, f))

        return size

    def get_files(self, path: str) -> typing.Tuple[str, float]:
        dirname = urllib_parse.urlparse(path).path
        dirname = self._to_namesafe(dirname)
        dirname = os.path.join(self._local_path, dirname)

        if os.path.isdir(dirname):
            size = self._download_size(dirname)
            if size == 0.0:
                size = self._download(path, dirname)

        else:
            os.mkdir(dirname)
            size = self._download(path, dirname)

        return dirname, size

    def _tar_files(self, dirs: typing.List[str], tarname: str):
        if os.path.isfile(tarname):
            return

        try:
            f = tarfile.open(tarname, 'w')
            for d in dirs:
                files = os.listdir(d)
                for fn in files:
                    f.add(os.path.join(d, fn))
            f.close()

        except Exception as e:
            os.remove(tarname)
            print('Error writing tar file %s. %s' % (tarname, e))

    def _make_name(self, dirs: typing.List[str]) -> typing.Tuple[str, str]:
        # name will be first and last entries
        d1 = self._to_namesafe(dirs[0].lstrip(self._local_path))
        d2 = self._to_namesafe(dirs[-1].lstrip(self._local_path))

        tarname = '%s-%s' % (d1, d2)
        uploadname = 'AFLOWLIB_%s' % tarname
        tarname = os.path.join(self._local_path, tarname + '.tar')

        return tarname, uploadname

    def _cleanup(self, ilist: typing.Union[str, typing.List[str]]):
        if isinstance(ilist, str):
            ilist = [ilist]
        for name in ilist:
            subprocess.Popen(['rm', '-rf', name])

    def _is_done_upload(self, uid: int) -> bool:
        res = self.client.uploads.get_upload(upload_id=uid).response().result
        Nproc = res.processed_calcs
        Ncalc = res.total_calcs

        if Nproc != Ncalc:
            return False
        return True

    def _get_status_upload(self, uploadname: str) -> typing.Tuple[str, int]:
        res = self.client.uploads.get_uploads(name=uploadname, state='all').response().result
        entries = res.results
        status = None
        upload_id = None

        for entry in entries:
            if entry['name'] == uploadname:
                status = 'uploaded'
                if entry['published']:
                    status = 'published'
                upload_id = entry['upload_id']
                break

        return status, upload_id

    def get_payload(self, uid: int) -> typing.Dict[str, typing.Any]:
        timenow = datetime.datetime.utcnow()
        if self.db_name == 'aflowlib':
            return dict(
                operation='publish',
                metadata=dict(
                    with_embargo=False,
                    comment='',
                    references=[
                        'http://www.sciencedirect.com/science/article/pii/S0927025614003322',
                        'http://aflowlib.org',
                        'http://www.sciencedirect.com/science/article/pii/S0927025612000687'],
                    coauthors=[
                        'f409d859-2639-4f82-b198-85e1c7c62f8b',
                        '580f7036-97b8-42b1-a9e6-815058bcac72',
                        '4c8d767d-335f-4ccd-9459-d0152b2026bc',
                        '7effd16a-a65c-4d95-b692-652488d94146',
                        'd59b3610-6335-4ad8-aca0-24a905de3a25',
                        '3df68fed-6ca0-4bf9-b2b1-d6be71e18f72',
                        'f7b540eb-a266-4379-9611-911ed8e3630e',
                        '7b3fe468-0011-4ba8-bd53-1ee65acda114',
                        '48e7a028-1a41-440d-9986-38540f5079c9',
                        'd63d07e6-ccc8-4eac-82af-d841627b6c53',
                        '9c308f66-eed1-4996-b33e-af78fb4944c7',
                        'd2621bc7-c45a-4d35-9dc1-5c05fa8326cb',
                        'ecba0e68-65ee-4b40-8fbf-a42714b1072b',
                        '81b96683-7170-49d7-8c4e-e9f34906b3ea'],
                    shared_with=[],
                    _upload_time=timenow,
                    _uploader='81b96683-7170-49d7-8c4e-e9f34906b3ea'))
        else:
            raise NotImplementedError('%s not yet supported.' % self.db_name)

    def publish(self, uids=None):
        print('Publishing')
        if uids is None:
            uids = self.uids
        for uid in uids:
            if self._is_done_upload(uid):
                payload = self.get_payload(uid)
                self.client.uploads.exec_upload_operation(upload_id=uid, payload=payload).response()

    def upload(self, file_path: str, name: str) -> int:
        res = self.client.uploads.upload(
            local_path=os.path.abspath(file_path), name=name).response().result
        return res.upload_id

    def download_proc(self, plist: typing.List[str]):
        size = 0.0
        max_zip_size = config.max_upload_size
        dirs = []
        done = []
        for i in range(len(plist)):
            d, s = self.get_files(self.update_list[plist[i]])
            if not self.do_upload:
                continue

            size += s
            dirs.append(d)
            done.append(plist[i])
            if size > max_zip_size or i == (len(plist) - 1):
                tarname, uploadname = self._make_name(dirs)
                status, uid = self._get_status_upload(uploadname)
                if status == 'published':
                    continue
                if status != 'uploaded':
                    self._tar_files(dirs, tarname)
                    uid = nomad_upload.upload_file(tarname, name=uploadname, offline=True)
                if self.do_publish:
                    self.publish([uid])
                self.uids.append(uid)
                if self.cleanup:
                    self._cleanup(dirs)
                    self._cleanup(tarname)
                size = 0.0
                dirs = []

    def download(self):
        print('Downloading from %s' % self.root_url)
        s = time.time()
        plist = [[] for i in range(self.parallel)]
        cur = 0
        for i in range(len(self.update_list)):
            if self.is_updated_list[i]:
                continue
            plist[cur % self.parallel].append(i)
            cur += 1

        if self.parallel > 1:
            procs = []
            for i in range(self.parallel):
                p = threading.Thread(target=self.download_proc, args=(plist[i],))
                procs.append(p)
            [p.start() for p in procs]
            [p.join() for p in procs]

        else:
            self.download_proc(plist[0])

        print('Time for download and upload (s)', time.time() - s)

    def prep_list(self):
        if self.outfile is not None and os.path.isfile(self.outfile):
            self.update_list = []
            self.is_updated_list = []
            for list_from_file in self._read_from_file(self.outfile):
                if isinstance(list_from_file, str):
                    self.update_list.append(list_from_file)
                    self.is_updated_list.append(False)
                else:
                    self.update_list.append(list_from_file[0])
                    self.is_updated_list.append(list_from_file[1])

        else:
            if self.parallel > 1:
                procs = []
                procs.append(threading.Thread(target=self.get_db_list))
                procs.append(threading.Thread(target=self.get_nomad_list))
                [p.start() for p in procs]
                [p.join() for p in procs]
            else:
                self.get_db_list()
                self.get_nomad_list()
            self.compare_lists()

    def update(self):
        self.prep_list()
        if self.do_download:
            self.download()


@client.command(
    help='Synchronizes the NOMAD database with the given external database.')
@click.argument('db_name', nargs=1, required=True)
@click.argument('root_url', nargs=1, required=True)
@click.option(
    '--outfile', default=None,
    help='File to read/write files missing in NOMAD database')
@click.option(
    '--nomadfile', default=None,
    help='File to read/write files in NOMAD database')
@click.option(
    '--dbfile', default=None,
    help='File to read/write files in given database')
@click.option(
    '--parallel', default=2,
    help='Number of processes to spawn to download/upload files')
@click.option(
    '--do-download', is_flag=True, default=False,
    help='Flag to automatically download downloaded files')
@click.option(
    '--do-upload', is_flag=True, default=False,
    help='Flag to automatically upload downloaded files')
@click.option(
    '--do-publish', is_flag=True, default=False,
    help='Flag to automatically publish upload')
@click.option(
    '--cleanup', is_flag=True, default=False,
    help='Flag to clean up downloaded files')
def synchdb(**kwargs):
    db = DbUpdater(**kwargs)
    db.update()
