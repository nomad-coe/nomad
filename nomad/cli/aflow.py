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

# TODO The code with changed to use the v1 API. This is not tested and most likely this
# code will fail.
# TODO The metadata should not be set via API, but added to the uploads as nomad.json.

from typing import List
import requests
import re
import subprocess
from urllib import parse as urllib_parse
import os
import tarfile
import threading
import time
import typing
import io
import re
import uuid
import json
import numpy as np
import ase
import bs4
import matid  # pylint: disable=import-error

from nomad import atomutils, config, client, processing as proc
from nomad.client import api, upload_file


class DbUpdater:
    '''
    Automatically synchronizes nomad it with a given database. It creates a list of paths
    to mainfiles in nomad and compares it with paths in the external database. The missing
    paths in nomad will then be downloaded from the external database and subsequently
    uploaded to nomad. The downloaded files are by default saved in '/nomad/fairdi/external'.
    '''

    def __init__(self, *args, auth: client.Auth, **kwargs):
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
        self.max_depth = None
        self.target_file = None
        self.uids: List[str] = []
        self._set(**kwargs)
        self.auth = auth

    def _set(self, **kwargs):
        for key, val in kwargs.items():
            key = key.lower()
            if hasattr(self, key):
                setattr(self, key, val)
            else:
                raise KeyError('Invalid key %s' % key)

        # create directory to save data
        subdir = ''
        if self.db_name.lower() == 'aflowlib':
            subdir = os.path.basename(self.root_url.strip('/'))

        dbpath = os.path.join(self.local_path, self.db_name)
        self._local_path = os.path.join(dbpath, subdir)
        if not os.path.isdir(dbpath):
            os.mkdir(dbpath)
        if not os.path.isdir(self._local_path):
            os.mkdir(self._local_path)

        # set target file based on database
        if self.db_name.lower() == 'aflowlib':
            if self.target_file is None:
                self.target_file = 'vasprun.xml.relax2.xz' if re.match(
                    r'.+/LIB\d+_LIB/?', self.root_url) else 'OUTCAR.static.xz'
        else:
            raise NotImplementedError('%s not yet supported.' % self.db_name)

        self._session = requests.Session()

    def _get_paths(self, root: str) -> typing.List[str]:
        response = self._session.get(root, verify=False)
        if not response.ok:
            response.raise_for_status()

        paths = []
        for url in re.findall('<a href="([^"]+)">', str(response.content)):
            if url in response.url:
                continue
            paths.append(os.path.join(response.url, url.lstrip('./')))
        return paths

    def _is_valid_file(self, path: str) -> bool:
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

    def get_db_list(self):
        if self.dbfile is not None and os.path.isfile(self.dbfile):
            self.db_files = self._read_from_file(self.dbfile)
        else:
            print('Generating list from %s' % self.root_url)
            self.db_files = []
            todo = self._get_paths(self.root_url)
            while len(todo) > 0:
                cur = todo.pop(0)
                depth = urllib_parse.urlparse(cur).path.rstrip('/').count('/')
                if self.target_file in cur and self.max_depth is None:
                    # it makes the assumption that the the mainfiles are located at the
                    # same level
                    self.max_depth = depth - 1

                if self.max_depth is not None and depth > self.max_depth:
                    pass

                elif self.max_depth is not None and depth == self.max_depth:
                    self.db_files.append(cur)

                else:
                    try:
                        add = self._get_paths(cur)
                    except Exception:
                        add = []
                    add = [path for path in add if self._is_valid_file(path)]
                    todo = add + todo

            if self.dbfile is not None:
                self._write_to_file(self.db_files, self.dbfile)

    def get_nomad_list(self):
        if self.nomadfile is not None and os.path.isfile(self.nomadfile):
            self.nomad_files = self._read_from_file(self.nomadfile)
        else:
            print('Generating NOMAD list')
            if self.db_name.lower() == 'aflowlib':
                servers = ['LIB%d_LIB' % n for n in range(1, 10)] + ['ICSD_WEB']
                paths = [s for s in servers if s in self.root_url]
                paths = paths if paths else servers
                # main_author: Stefano Curtarolo
                query = dict(
                    main_author_id='81b96683-7170-49d7-8c4e-e9f34906b3ea',
                    paths=paths)

            self.nomad_files = []
            page_after_value = None
            while True:
                response = api.post('entries/query', data=json.dumps({
                    'owner': 'all',
                    'query': query,
                    'aggregations': {
                        'mainfiles': {
                            'terms': {
                                'quantity': 'mainfile',
                                'pagination': {
                                    'page_size': 1000,
                                    'page_after_value': page_after_value
                                }
                            }
                        }
                    }
                }))
                assert response.status_code == 200
                aggregation = response.json()['aggregations']['mainfiles']['terms']
                page_after_value = aggregation['pagination']['next_page_after_value']
                for bucket in aggregation['data']:
                    self.nomad_files.append(bucket['value'])

                if len(aggregation['data']) < 1000 or page_after_value is None:
                    break

            if self.nomadfile is not None:
                self._write_to_file(self.nomad_files, self.nomadfile)

    def compare_lists(self):
        '''
        Identify the difference between the nomad list and db list
        '''
        def reduce_list(ilist: typing.List[str]):
            olist = []
            for e in ilist:
                p = urllib_parse.urlparse(e).path.strip('/')
                olist.append(os.path.join(*p.split('/')[1:self.max_depth]))
            olist = list(set(olist))
            olist.sort()
            return olist

        print('Identifying differences')
        self.max_depth = 10 if self.max_depth is None else self.max_depth
        db = reduce_list(self.db_files)
        nomad = reduce_list(self.nomad_files)
        ns = set(nomad)
        self.update_list = [i for i in db if i not in ns]

        ds = set(db)
        in_nomad = [i for i in nomad if i not in ds]
        if len(in_nomad) > 0:
            fn = 'in_nomad.txt'
            print('Warning: Some NOMAD entries not found in db.')
            print('See %s for list.' % fn)
            self._write_to_file(in_nomad, fn)

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

    def _get_files(self, path: str) -> typing.Tuple[str, float]:
        def is_dir(path: str) -> bool:
            path = path.strip()
            if path[-1] == '/' and self.root_url in path:
                return True
            return False

        def download(path: str, iodir: str) -> float:
            files = self._get_paths(path)
            files = [f for f in files if self._is_valid_file(f)]
            size = 0.0
            for f in files:
                if is_dir(f):
                    _, s = self._get_files(f)
                    size += s
                    continue

                res = self._session.get(f, stream=True)
                fn = res.url.split('/')[-1]
                fn = os.path.join(iodir, fn)
                try:
                    with open(fn, 'wb') as fb:
                        fb.write(res.content)
                except Exception:
                    continue

                size += os.path.getsize(fn)
            return size

        def get_download_size(dirname: str) -> float:
            complete = False
            files = os.listdir(dirname)

            if self.db_name.lower() == 'aflowlib':
                complete = self.target_file in files
            else:
                complete = True

            size = 0.0
            if not complete:
                self._cleanup([os.path.join(dirname, f) for f in files])
            else:
                for f in files:
                    size += os.path.getsize(os.path.join(dirname, f))

            return size

        dirname = urllib_parse.urlparse(path).path
        dirname = self._to_namesafe(dirname)
        dirname = os.path.join(self._local_path, dirname)

        if os.path.isdir(dirname):
            size = get_download_size(dirname)
            if size == 0.0:
                size = download(path, dirname)

        else:
            try:
                os.mkdir(dirname)
                size = download(path, dirname)
            except Exception:
                size = 0.0

        return dirname, size

    def _make_name(self, dirs: typing.List[str]) -> typing.Tuple[str, str]:
        # name will be first and last entries
        d1 = self._to_namesafe(dirs[0].lstrip(self._local_path))
        d2 = self._to_namesafe(dirs[-1].lstrip(self._local_path))

        tarname = '%s-%s' % (d1, d2)
        uploadname = '%s_%s' % (self.db_name.upper(), tarname)
        tarname = os.path.join(self._local_path, '%s.tar' % tarname)

        return tarname, uploadname

    def _cleanup(self, ilist: typing.Union[str, typing.List[str]]):
        if isinstance(ilist, str):
            ilist = [ilist]
        for name in ilist:
            subprocess.Popen(['rm', '-rf', name])

    # TODO This should not be set via API, but added as nomad.json to the upload!
    # def get_payload(self, uid: int) -> typing.Dict[str, typing.Any]:
    #     if self.db_name == 'aflowlib':
    #         return dict(
    #             operation='publish',
    #             metadata=dict(
    #                 with_embargo=False,
    #                 comment='',
    #                 references=[
    #                     'http://www.sciencedirect.com/science/article/pii/S0927025614003322',
    #                     'http://aflowlib.org',
    #                     'http://www.sciencedirect.com/science/article/pii/S0927025612000687'],
    #                 coauthors=[
    #                     'f409d859-2639-4f82-b198-85e1c7c62f8b',
    #                     '580f7036-97b8-42b1-a9e6-815058bcac72',
    #                     '4c8d767d-335f-4ccd-9459-d0152b2026bc',
    #                     '7effd16a-a65c-4d95-b692-652488d94146',
    #                     'd59b3610-6335-4ad8-aca0-24a905de3a25',
    #                     '3df68fed-6ca0-4bf9-b2b1-d6be71e18f72',
    #                     'f7b540eb-a266-4379-9611-911ed8e3630e',
    #                     '7b3fe468-0011-4ba8-bd53-1ee65acda114',
    #                     '48e7a028-1a41-440d-9986-38540f5079c9',
    #                     'd63d07e6-ccc8-4eac-82af-d841627b6c53',
    #                     '9c308f66-eed1-4996-b33e-af78fb4944c7',
    #                     'd2621bc7-c45a-4d35-9dc1-5c05fa8326cb',
    #                     'ecba0e68-65ee-4b40-8fbf-a42714b1072b',
    #                     '81b96683-7170-49d7-8c4e-e9f34906b3ea'],
    #                 shared_with=[]))
    #     else:
    #         raise NotImplementedError('%s not yet supported.' % self.db_name)

    def publish(self, uids=None):
        def is_done_upload(uid: int) -> bool:
            response = api.get(f'uploads/{uid}', auth=self.auth)
            assert response.status_code == 200
            return response.json()['data']['process_status'] in proc.ProcessStatus.STATUSES_NOT_PROCESSING

        print('Publishing')
        uids = self.uids if uids is None else uids
        for uid in uids:
            if is_done_upload(uid):
                response = api.post(f'uploads/{uid}/action/publish', auth=self.auth)
                assert response.status_code == 200

    def upload(self, file_path: str, upload_name: str) -> int:
        uid = upload_file(os.path.abspath(file_path), self.auth, local_path=True, upload_name=upload_name)
        assert uid is not None
        return uid

    def _download_proc(self, plist: typing.List[str]):
        def tar_files(dirs: typing.List[str], tarname: str):
            if os.path.isfile(tarname):
                return

            try:
                with tarfile.open(tarname, 'w') as f:
                    for d in dirs:
                        files = os.listdir(d)
                        for fn in files:
                            f.add(os.path.join(d, fn))

            except Exception as e:
                os.remove(tarname)
                print('Error writing tar file %s. %s' % (tarname, e))

        def get_status_upload(uploadname: str) -> typing.Tuple[str, str]:
            response = api.get(f'uploads', params=dict(name=uploadname), auth=self.auth)
            assert response.status_code == 200
            response_json = response.json()
            assert len(response_json['data']) <= 1

            if len(response_json['data']) == 0:
                return None, None

            upload = response_json['data'][0]
            if upload['published']:
                return 'published', upload['upload_id']

            return 'uploaded', upload['upload_id']

        size = 0.0
        max_zip_size = config.process.max_upload_size
        dirs = []
        for i in range(len(plist)):
            d, s = self._get_files(self.update_list[plist[i]])
            if not self.do_upload:
                continue

            size += s
            dirs.append(d)
            if size > max_zip_size or i == (len(plist) - 1):
                tarname, uploadname = self._make_name(dirs)
                status, uid = get_status_upload(uploadname)
                if status == 'published':
                    continue
                if status != 'uploaded':
                    tar_files(dirs, tarname)
                    uid = upload_file(tarname, auth=self.auth, upload_name=uploadname, local_path=True)
                    assert uid is not None
                if self.do_publish:
                    self.publish([uid])
                self.uids.append(uid)
                if self.cleanup:
                    self._cleanup(dirs)
                    self._cleanup(tarname)
                size = 0.0
                dirs = []

    def download(self):
        '''
        Download files from database.
        '''
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
                p = threading.Thread(target=self._download_proc, args=(plist[i],))
                procs.append(p)
            [p.start() for p in procs]
            [p.join() for p in procs]

        else:
            self._download_proc(plist[0])

        print('Time for download and upload (s)', time.time() - s)

    def get_list_to_download(self):
        '''
        Generate lists of files from database and from nomad and returns the difference.
        '''
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
        self.get_list_to_download()
        if self.do_download:
            self.download()


def write_prototype_data_file(aflow_prototypes: dict, filepath) -> None:
    '''Writes the prototype data file in a compressed format to a python
    module.

    Args:
        aflow_prototypes
    '''
    class NoIndent(object):
        def __init__(self, value):
            self.value = value

    class NoIndentEncoder(json.JSONEncoder):
        '''A custom JSON encoder that can pretty-print objects wrapped in the
        NoIndent class.
        '''
        def __init__(self, *args, **kwargs):
            super(NoIndentEncoder, self).__init__(*args, **kwargs)
            self.kwargs = dict(kwargs)
            del self.kwargs['indent']
            self._replacement_map = {}

        def default(self, o):  # pylint: disable=E0202
            if isinstance(o, NoIndent):
                key = uuid.uuid4().hex
                self._replacement_map[key] = json.dumps(o.value, **self.kwargs)
                return "@@%s@@" % (key,)
            else:
                return super(NoIndentEncoder, self).default(o)

        def encode(self, o):
            result = super(NoIndentEncoder, self).encode(o)
            for k, v in self._replacement_map.items():
                result = result.replace('"@@%s@@"' % (k,), v)
            return result

    prototype_dict = aflow_prototypes["prototypes_by_spacegroup"]
    for prototypes in prototype_dict.values():
        for prototype in prototypes:
            # Save the information back in a prettified form
            prototype["atom_positions"] = NoIndent(prototype["atom_positions"])
            prototype["atom_labels"] = NoIndent(prototype["atom_labels"])
            prototype["lattice_vectors"] = NoIndent(prototype["lattice_vectors"])
            try:
                prototype["normalized_wyckoff_matid"] = NoIndent(prototype["normalized_wyckoff_matid"])
            except KeyError:
                pass

    # Save the updated data
    with io.open(filepath, "w", encoding="utf8") as f:
        json_dump = json.dumps(aflow_prototypes, ensure_ascii=False, indent=4, sort_keys=True, cls=NoIndentEncoder)
        json_dump = re.sub(r"\"(-?\d+(?:[\.,]\d+)?)\"", r'\1', json_dump)  # Removes quotes around numbers
        f.write("# -*- coding: utf-8 -*-\naflow_prototypes = {}\n".format(json_dump))


def update_prototypes(ctx, filepath, matches_only):

    if matches_only:
        from nomad.aflow_prototypes import aflow_prototypes
    else:
        # The basic AFLOW prototype data is available in a Javascript file. Here we
        # retrieve it and read only the prototype list from it.
        prototypes_file_url = 'http://aflowlib.org/CrystalDatabase/js/table_sort.js'
        r = requests.get(prototypes_file_url, allow_redirects=True)
        datastring = r.content.decode("utf-8")
        datastring = datastring.split('];')[0]
        datastring = datastring.split('= [')[1]
        data = json.loads('[' + datastring + ']')

        newdictarray = []
        n_prototypes = 0
        n_missing = 0
        for protodict in data:
            n_prototypes += 1
            newdict = {}

            # Make prototype plaintext
            prototype = bs4.BeautifulSoup(protodict["Prototype"], "html5lib").getText()

            # Add to new dictionary
            newdict['Notes'] = protodict['Notes']
            newdict['Prototype'] = prototype
            newdict['Space Group Symbol'] = protodict['Space Group Symbol']
            newdict['Space Group Number'] = protodict['Space Group Number']
            newdict['Pearsons Symbol'] = protodict['Pearson Symbol']
            newdict['Strukturbericht Designation'] = protodict['Strukturbericht Designation']
            newdict['aflow_prototype_id'] = protodict['AFLOW Prototype']
            newdict['aflow_prototype_url'] = 'http://www.aflowlib.org/CrystalDatabase/' + protodict['href'][2:]

            # Download cif or poscar if possible make ASE ase.Atoms object if possible
            # to obtain labels, positions, cell
            cifurl = 'http://www.aflowlib.org/CrystalDatabase/CIF/' + protodict['href'][2:-5] + '.cif'
            r = requests.get(cifurl, allow_redirects=True)
            cif_str = r.content.decode("utf-8")
            cif_file = io.StringIO()
            cif_file.write(cif_str)
            cif_file.seek(0)
            try:
                atoms = ase.io.read(cif_file, format='cif')
            except Exception:
                print("Error in getting prototype structure from CIF: {}", format(cifurl))
                # Then try to get structure from POSCAR
                try:
                    poscarurl = 'http://www.aflowlib.org/CrystalDatabase/POSCAR/' + protodict['href'][2:-5] + '.poscar'
                    r = requests.get(poscarurl, allow_redirects=True)
                    poscar_str = r.content.decode("utf-8")
                    poscar_file = io.StringIO()
                    poscar_file.write(poscar_str)
                    poscar_file.seek(0)
                    atoms = ase.io.read(poscar_file, format='vasp')
                except Exception:
                    print("Error in getting prototype structure from POSCAR: {}".format(poscarurl))
                    print("Could not read prototype structure from CIF or POSCAR file for prototype: {}, {}, ".format(prototype, newdict['aflow_prototype_url']))
                    n_missing += 1
                    continue

            atom_positions = atoms.get_positions()
            atom_labels = atoms.get_chemical_symbols()
            cell = atoms.get_cell()

            newdict['lattice_vectors'] = cell.tolist()
            newdict['atom_positions'] = atom_positions.tolist()
            newdict['atom_labels'] = atom_labels
            newdictarray.append(newdict)

            print("Processed: {}".format(len(newdictarray)))

        # Sort prototype dictionaries by spacegroup and make dictionary
        structure_types_by_spacegroup = {}
        for i_sg in range(1, 231):
            protos_sg = []
            for newdict in newdictarray:
                if newdict['Space Group Number'] == i_sg:
                    protos_sg.append(newdict)
            structure_types_by_spacegroup[i_sg] = protos_sg

        # Wrap in a dictionary that can hold other data, e.g. the symmemtry tolerance parameter.
        aflow_prototypes = {
            "prototypes_by_spacegroup": structure_types_by_spacegroup
        }
        print(
            "Extracted latest AFLOW prototypes online. Total number of "
            "successfully fetched prototypes: {}, missing: {}"
            .format(n_prototypes, n_missing)
        )

    # Update matches
    n_prototypes = 0
    n_failed = 0
    n_unmatched = 0
    prototype_dict = aflow_prototypes["prototypes_by_spacegroup"]

    for aflow_spg_number, prototypes in prototype_dict.items():
        n_prototypes += len(prototypes)
        for prototype in prototypes:

            # Read prototype structure
            pos = np.array(prototype["atom_positions"])
            labels = prototype["atom_labels"]
            cell = np.array(prototype["lattice_vectors"])
            atoms = ase.Atoms(
                symbols=labels,
                positions=pos,
                cell=cell,
                pbc=True
            )

            # Try to first see if the space group can be matched with the one in AFLOW
            try:
                symm = matid.SymmetryAnalyzer(atoms, config.normalize.prototype_symmetry_tolerance)
                spg_number = symm.get_space_group_number()
                wyckoff_matid = symm.get_wyckoff_letters_conventional()
                norm_system = symm.get_conventional_system()
            except Exception:
                n_failed += 1
            else:
                # If the space group is matched, add the MatID normalized Wyckoff
                # letters to the data.
                if spg_number == aflow_spg_number:
                    atomic_numbers = norm_system.get_atomic_numbers()
                    normalized_wyckoff_matid = atomutils.get_normalized_wyckoff(atomic_numbers, wyckoff_matid)
                    prototype["normalized_wyckoff_matid"] = normalized_wyckoff_matid
                else:
                    n_unmatched += 1
    print(
        "Updated matches in AFLOW prototype library. Total number of "
        "prototypes: {}, unmatched: {}, failed: {}"
        .format(n_prototypes, n_unmatched, n_failed)
    )

    # Write data file to the specified path
    write_prototype_data_file(aflow_prototypes, filepath)
