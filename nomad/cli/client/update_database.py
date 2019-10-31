import requests
import re
import subprocess
from urllib.parse import urlparse
import os
from .client import client
from nomad import config
from nomad.cli.client import upload as nomad_upload
import datetime
import click
import tarfile
import threading
import time


class DbUpdater:
    def __init__(self, *args, **kwargs):
        self.db_name = args[0].strip()
        self.root_url = args[1].strip()
        self.user = 'admin'
        self.password = 'password'
        self.local_path = '/nomad/fairdi/external'
        self.nomad_server = 'http://labdev-nomad.esc.rzg.mpg.de/fairdi/nomad/testing/api/uploads/'
        self.max_zip_size = 32000000000
        self.do_download = False
        self.do_upload = False
        self.do_publish = False
        self.nproc = 2
        self.uids = []
        self._set(**kwargs)
        self._set_db()
        self._set_local_path()
        self.configure_client()

    def _set(self, **kwargs):
        for key, val in kwargs.items():
            if key.lower() == 'user':
                self.user = val
            elif key.lower() == 'password':
                self.password = val
            elif key.lower() == 'local_path':
                self.local_path = val
            elif key.lower() == 'nomad_server':
                self.nomad_server = val
            elif key.lower() == 'outfile':
                self.outfile = val
            elif key.lower() == 'nomadfile':
                self.nomadfile = val
            elif key.lower() == 'dbfile':
                self.dbfile = val
            elif key.lower() == 'nproc':
                self.nproc = val
            elif key.lower() == 'download':
                self.do_download = val
            elif key.lower() == 'upload':
                self.do_upload = val
            elif key.lower() == 'publish':
                self.do_publish = val
            else:
                raise KeyError('Invalid key %s' % key)

    def configure_client(self, **kwargs):
        from nomad.cli.client import create_client
        i = self.nomad_server.index('api')
        with open('nomad.yaml', 'w') as f:
            f.write('client:\n')
            f.write('    user: %s\n' % self.user)
            f.write('    password: %s\n' % self.password)
            f.write('    url: %s\n' % self.nomad_server[:i + 3])
        config.client.url = self.nomad_server[:i + 3]
        config.client.user = self.user
        config.client.password = self.password
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

    def _open_page(self, path):
        with requests.Session() as session:
            response = session.get(path, verify=False)
        if not response.ok:
            return response.raise_for_status()
        return response

    def get_paths(self, root):
        response = self._open_page(root)
        paths = []
        for url in re.findall('<a href="([^"]+)">', str(response.content)):
            if url in response.url:
                continue
            paths.append(os.path.join(response.url, url.lstrip('./')))
        return paths

    def _filter_files(self, files):
        new = []
        for f in files:
            if self._is_file(f):
                new.append(f)
        return new

    def _is_file(self, path):
        if 'vasprun.xml' in path:
            return True
        return False

    def _is_dir(self, path):
        path = path.strip()
        if path[-1] == '/' and self.root_url in path:
            return True
        return False

    def _depth(self, path):
        return urlparse(path).path.rstrip('/').count('/')

    def _rules_ok(self, path):
        ok = False
        if self._is_file(path):
            ok = True
        if self._depth(path) >= self.max_depth:
            ok = True
        if '?' in path:
            ok = False
        return ok

    def _to_string(self, path):
        return path.strip('/').replace('/', '_').replace(':', '_')

    def _read_from_file(self, filename):
        print('Reading from file %s' % filename)
        with open(filename) as f:
            data = f.readlines()
        data = [s.strip().split() for s in data]
        for i in range(len(data)):
            for j in range(len(data[i])):
                data[i][j] = data[i][j].strip()
                if data[i][j] in ['T', 'True', '1', 'Y']:
                    data[i][j] = True
                elif data[i][j] in ['F', 'False', '0', 'N']:
                    data[i][j] = False
            if len(data[i]) == 1:
                data[i] = data[i][0]
        return data

    def _write_to_file(self, data, filename):
        with open(filename, 'w') as f:
            for i in range(len(data)):
                if isinstance(data[i], str):
                    f.write('%s\n' % data[i])
                else:
                    for j in range(len(data[i])):
                        e = data[i][j]
                        if e is True:
                            e = 'T'
                        elif e is False:
                            e = 'F'
                        else:
                            e = str(e)
                        f.write('%s ' % e)
                    f.write('\n')

    def _gen_db_list(self):
        print('Generating list from %s' % self.root_url)
        self.db_files = []
        todo = self.get_paths(self.root_url)
        while len(todo) > 0:
            cur = todo[-1]
            if self._rules_ok(cur):
                self.db_files.append(cur)
            elif self._is_dir(cur):
                add = self.get_paths(cur)
                todo = add + todo
            todo.pop(-1)
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
            servers = ['LIB1_LIB', 'LIB2_LIB', 'LIB3_LIB']
            paths = []
            for s in servers:
                if s in self.root_url:
                    paths.append(s)
            if len(paths) == 0:
                paths = servers
            kwargs = dict(authors=['Curtarolo, Stefano'], paths=paths, scroll=True)

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
            for i in range(len(self.nomad_files)):
                print(self.nomad_files[i])
            self._write_to_file(self.nomad_files, self.nomadfile)

    def get_nomad_list(self):
        if self.nomadfile is not None and os.path.isfile(self.nomadfile):
            self.nomad_files = self._read_from_file(self.nomadfile)
        else:
            self._gen_nomad_list()

    def _reduce_list(self, ilist):
        olist = []
        for e in ilist:
            p = urlparse(e).path.strip('/')
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
        u = urlparse(self.root_url)
        up = u.path.strip('/').split('/')[0]
        root = '%s://%s/%s' % (u.scheme, u.netloc, up)
        self.update_list = [os.path.join(root, e) for e in self.update_list]
        self.is_updated_list = [False] * len(self.update_list)
        print('Found %d entries to be added in NOMAD' % len(self.update_list))
        if self.outfile is not None:
            data = [self.update_list[i] for i in range(len(self.update_list))]
            self._write_to_file(data, self.outfile)

    def _download(self, path, iodir):
        files = self.get_paths(path)
        files = [f for f in files if self._rules_ok(f)]
        size = 0.0
        with requests.Session() as session:
            for f in files:
                res = session.get(f, stream=True)
                fn = res.url.split('/')[-1]
                fn = os.path.join(iodir, fn)
                with open(fn, 'wb') as fb:
                    fb.write(res.content)
                size += os.path.getsize(fn)
        return size

    def get_files(self, path):
        dirname = urlparse(path).path
        dirname = self._to_string(dirname)
        dirname = os.path.join(self._local_path, dirname)
        if os.path.isdir(dirname):
            size = 0.0
            files = os.listdir(dirname)
            complete = False
            for f in files:
                if 'vasprun' in f:
                    complete = True
                size += os.path.getsize(os.path.join(dirname, f))
            if not complete:
                self.cleanup([os.path.join(dirname, f) for f in files])
                size = self._download(path, dirname)
        else:
            os.mkdir(dirname)
            size = self._download(path, dirname)
        return dirname, size

    def _tar_files(self, dirs, tarname):
        if os.path.isfile(tarname):
            return
        with tarfile.open(tarname, 'w') as f:
            for d in dirs:
                files = os.listdir(d)
                for fn in files:
                    f.add(os.path.join(d, fn))

    def cleanup(self, ilist):
        if isinstance(ilist, str):
            ilist = [ilist]
        for name in ilist:
            subprocess.Popen(['rm', '-rf', name])

    def _is_done_upload(self, uid):
        res = self.client.uploads.get_upload(upload_id=uid).response().result
        Nproc = res.processed_calcs
        Ncalc = res.total_calcs
        if Nproc != Ncalc:
            return False
        return True

    def get_payload(self, uid):
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
                        148, 149, 150, 151, 152, 146, 145, 138, 137,
                        136, 135, 134, 147, 125],
                    shared_with=[],
                    _upload_time=timenow,
                    _uploader=125))

    def publish(self):
        print('Publishing')
        for uid in self.uids:
            if self._is_done_upload(uid):
                payload = self.get_payload(uid)
                self.client.uploads.exec_upload_operation(upload_id=uid, payload=payload)

    def register_upload(self, donelist, plist, pn):
        data = []
        for i in plist:
            if i in donelist:
                self.is_updated_list[i] = True
            else:
                data.append(self.update_list[i])
        self._write_to_file(data, self.outfile + '_%d' % pn)

    def aggregate_procs(self):
        data = []
        for i in range(self.nproc):
            data += self._read_from_file(self.outfile + '_%d' % i)
        self._write_to_file(data, self.outfile + '_updated')

    def download_proc(self, plist, pn):
        size = 0.0
        dirs = []
        done = []
        for i in range(len(plist)):
            d, s = self.get_files(self.update_list[plist[i]])
            if not self.do_upload:
                continue
            size += s
            dirs.append(d)
            done.append(plist[i])
            if size > self.max_zip_size or i == (len(plist) - 1):
                tstamp = datetime.datetime.now().strftime('_%y%m%d%H%M%S%f')
                tname = self._to_string(dirs[0].lstrip(self._local_path))
                tname = os.path.join(self._local_path, tname + '%s.tar' % tstamp)
                self._tar_files(dirs, tname)
                uid = nomad_upload.upload_file(tname, name='AFLOWLIB_%s' % tstamp)
                self.uids.append(uid)
                self.register_upload(done, plist, pn)
                self.cleanup(dirs)
                self.cleanup(tname)
                size = 0.0
                dirs = []
                done = []
                return

    def download(self):
        print('Downloading from %s' % self.root_url)
        s = time.time()
        plist = [[] for i in range(self.nproc)]
        cur = 0
        for i in range(len(self.update_list)):
            if self.is_updated_list[i]:
                continue
            plist[cur % self.nproc].append(i)
            cur += 1
        procs = []
        for i in range(self.nproc):
            p = threading.Thread(target=self.download_proc, args=(plist[i], i,))
            procs.append(p)
        [p.start() for p in procs]
        [p.join() for p in procs]
        self.aggregate_procs()
        print('Time for download and upload (s)', time.time() - s)

    def prep_list(self):
        if self.outfile is not None and os.path.isfile(self.outfile):
            uplist = self._read_from_file(self.outfile)
            self.update_list = []
            self.is_updated_list = []
            for l in uplist:
                if isinstance(l, str):
                    self.update_list.append(l)
                    self.is_updated_list.append(False)
                else:
                    self.update_list.append(l[0])
                    self.is_updated_list.append(l[1])
        else:
            procs = []
            procs.append(threading.Thread(target=self.get_db_list))
            procs.append(threading.Thread(target=self.get_nomad_list))
            [p.start() for p in procs]
            [p.join() for p in procs]
            self.compare_lists()

    def update(self):
        self.prep_list()
        if self.do_download:
            self.download()
        if self.do_publish:
            self.publish()

@client.command(
    help='Synchronizes the NOMAD database with the given external database.')
@click.argument('database', nargs=1, required=True)
@click.argument('url', nargs=1, required=True)
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
    '--nproc', default=2,
    help='Number of processes to spawn to download/upload files')
@click.option(
    '--download', is_flag=True, default=False,
    help='Flag to automatically download downloaded files')
@click.option(
    '--upload', is_flag=True, default=False,
    help='Flag to automatically upload downloaded files')
@click.option(
    '--publish', is_flag=True, default=False,
    help='Flag to automatically publish upload')
def synchdb(database, url, outfile=None, nomadfile=None, dbfile=None, nproc=2, download=False, upload=False, publish=False):
    db = DbUpdater(database, url, outfile=outfile, nomadfile=nomadfile, dbfile=dbfile, nproc=nproc, download=download, upload=upload, publish=publish)
    db.update()
