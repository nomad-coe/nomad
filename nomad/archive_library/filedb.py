# Copyright 2019  Alvin Noe Ladines, Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Module for storage of archive data using the msgpack module.

In module ``ArchiveFileDB, the data are fragmented and saved as a list
in a msgpack file.
A component can be retrieved by giving an offset to the msgpack file to
be unpacked. The database can then be queried by giving a schema similar
to the archive data.

To build the database,

.. code-block:: python
    db = ArchiveFileDB("db.msg", mode='w', max_lfragment=3)
    db.add_data(["archive1.json", "archive2.json"])
    db.close()

To query the database,

.. code-block:: python
    db = ArchiveFileDB("db.msg", mode='r')
    db.query({'idX':{'sectionX':{'propertyX':None}}})
    db.close()
"""

import msgpack
import json
import os
from io import StringIO, BufferedWriter, BufferedReader, BytesIO


class JSONdata:
    """
    Provides a graphQL-style query for a given json data and query schema.
    Arguments:
        data: The json data to be queried
    """
    def __init__(self, data):
        self.data = data

    def _get_data(self, key, val):
        if isinstance(key, str):
            key = key.split('/')

        index = None
        bn0 = key[0].find('[')
        if bn0 > 0:
            bn1 = key[0].find(']')
            index = key[0][bn0 + 1:bn1]
            key[0] = key[0].replace(key[0][bn0:bn1 + 1], '')
            if ':' in index:
                v = val[key[0]]
                lo, hi = index.split(':')
                lo = int(lo) if lo else 0
                hi = int(hi) if hi else 0
                lo = len(v) + lo if lo < 0 else lo
                hi = len(v) + hi if hi <= 0 else hi
                if hi > len(v):
                    hi = len(v)
                if lo > hi or lo < 0 or hi < 0:
                    return
                index = list(range(lo, hi))
            else:
                index = int(index)

        if not key[0] in val:
            return

        v = val[key[0]]
        if isinstance(index, int):
            try:
                v = v[index]
            except (IndexError, KeyError, TypeError):
                return
        elif isinstance(index, list):
            try:
                v = [v[i] for i in index]
            except (IndexError, KeyError, TypeError):
                return
        if len(key) > 1:
            if isinstance(v, list):
                return [self._get_data(key[1:], vi) for vi in v]
            else:
                return self._get_data(key[1:], v)
        else:
            return v

    def get_data(self, entry, root=None):
        """
        Recursively searches the json data to fill in the null values
        in the given schema.
        Arguments:
            entry: a dict with a structure similar to the json data but
                with null value to be filled with the desired value.
        """
        data = {}
        if entry is None:
            v = self._get_data(root, self.data)
            return v
        else:
            for key, val in entry.items():
                if root is not None:
                    k = os.path.join(root, key)
                else:
                    k = key
                data[key] = self.get_data(val, k)
            return data


class ArchiveFileDB:
    """
    An interface to the messagepack module to provide an searchable
    container of archive data.
    Arguments:
        fileio: can be a string or file object to read/write the msgpack file
        mode: r/w to indicate read or write the msgpack file
        max_lfragment: the maximum level for which the archive data will
            be fragmented for more efficient unpacking of msgpack components
    """
    def __init__(self, fileio, mode='r', max_lfragment=None):
        self._filename = None
        self._fileobj = None
        if isinstance(fileio, str):
            self._filename = fileio
            self._mode = mode
        elif isinstance(fileio, BufferedReader):
            self._fileobj = fileio
            self._mode = 'rb'
        elif isinstance(fileio, BufferedWriter):
            self._fileobj = fileio
            self._mode = 'wb'
        elif isinstance(fileio, BytesIO):
            self._fileobj = fileio
            self._mode = mode
        else:
            raise TypeError
        self._max_lfragment = max_lfragment
        if 'w' in self._mode and self._max_lfragment is None:
            self._max_lfragment = 2
        if '+' in self._mode:
            self._max_lfragment = None
        self._sep = 'MSG_ENTRY'
        self._ids = None
        self._data = {}

    @property
    def max_lfragment(self):
        if self._max_lfragment is None:
            orig_mode = self.mode
            self.mode = 'rb'
            self._max_lfragment = self.get_docs('MAX_LFRAGMENT')
            self.mode = orig_mode
        return self._max_lfragment

    def _fragment_json(self, data, key='', cur_lfragment=0):
        if cur_lfragment >= self.max_lfragment or not isinstance(data, dict):
            return [dict(path=key, data={os.path.basename(key): data})]
        res = []
        cur_lfragment += 1
        main = dict(path=key, data=[])
        for k, v in data.items():
            p = os.path.join(key, k)
            res += self._fragment_json(v, p, cur_lfragment)
            main['data'].append(p)
        res += [main]
        return res

    def write(self, abspath, relpath):
        """
        Mimic the zipfile function to write files to database.
        Arguments:
            abspath: The absolute path to the file to be read
            relpath: For compatibility with zipfile
        """
        self.add_data(abspath)

    def close(self, save=True):
        """
        Mimic the zipfile function to close the msgpack file.
        Will trigger the creation of the database when in write mode.
        Arguments:
            save: If True will add the current data in memory to database
        """
        if 'w' in self._mode:
            self.create_db()
        if self._fileobj:
            self._fileobj.close()
            self._fileobj = None

    def save(self):
        """
        Commits current data in memory to database
        """
        self.create_db()

    def add_data(self, data):
        """
        Add data to the msgpack database.
        Arguments:
            data: Can be a filename or dictionary or list of both
        """
        if isinstance(data, str):
            key = os.path.basename(data)
            if data.endswith('json'):
                key = key.split('.')[0]
                self._data[key] = json.load(open(data))
            else:
                key = key.replace('.', '_')
                self._data[key] = open(data).read()
        elif isinstance(data, dict):
            key = list(data.keys())
            assert len(key) == 1
            self._data[key[-1]] = data[key[-1]]
        elif isinstance(data, list):
            for i in range(len(data)):
                self.add_data(data[i])
        else:
            raise NotImplementedError

    # TODO
    # def edit_data(self, data):
    #     if self._data:
    #         print("Memory not empty. Commit existing data first")
    #         return
    #     self.add_data(data)
    #     entries = self._fragment_json(self._data)
    #     data_to_write = self._load_data()
    #     raise NotImplementedError

    def _load_data(self):
        orig_mode = self.mode
        self.mode = 'rb'
        self.fileobj.seek(0)
        data_loaded = msgpack.load(self.fileobj)
        self.mode = orig_mode
        return data_loaded

    def create_db(self):
        """
        Creates the database and writes it to the msgpack file.
        The database consists of the list of the fragmented data
        and the list of footers such as the ids of the data.
        """
        data = self._data
        entries = self._fragment_json(data)
        # save the data in a separate list everytime function is called
        # last list reserved for IDs
        if '+' in self._mode:
            data_to_write = self._load_data()
            data_to_write.insert(-1, [])
        else:
            data_to_write = [[], []]

        # get last entry in database
        last_index = 0
        for i in range(len(data_to_write) - 1):
            last_index += len(data_to_write[i]) / 2
        last_index = int(last_index)

        for i in range(len(entries)):
            sep = '%s_%d' % (self._sep, i + last_index)
            data_to_write[-2].append(sep)
            data_to_write[-2].append(entries[i]['data'])
        data_str = msgpack.dumps(data_to_write)
        if not data_to_write[-1]:
            head = {}
        else:
            sep = '%s_IDS' % (self._sep)
            hi = data_to_write[-1].index(sep.encode()) + 1
            head = data_to_write[-1][hi]
            data_to_write[-1] = []
        sep = '%s_MAX_LFRAGMENT' % (self._sep)
        data_to_write[-1].append(sep)
        data_to_write[-1].append(self.max_lfragment)
        # add pointers to entries
        start_index = 0
        for i in range(len(entries)):
            sep = '%s_%d' % (self._sep, i + last_index)
            index = data_str.index(sep.encode(), start_index) + len(sep)
            head[entries[i]['path']] = index
            start_index = index
        sep = '%s_IDS' % (self._sep)
        data_to_write[-1].append(sep)
        data_to_write[-1].append(head)
        data_str = msgpack.dumps(data_to_write)
        self.fileobj.write(data_str)
        self._data = {}

    def _reduce_to_section(self, entry, cur_lfragment=0):
        if entry is None:
            return
        cur_lfragment += 1
        if cur_lfragment > self.max_lfragment:
            return
        new_dict = {}
        for key, val in entry.items():
            if '[' in key and ']' in key:
                key = key.split('[')[0]
            v = self._reduce_to_section(val, cur_lfragment)
            new_dict[key] = v
        return new_dict

    @staticmethod
    def to_list_path_str(entries, root='', paths=[]):
        if entries is None:
            return
        if len(paths) > 0:
            paths.remove(root)
        for key, val in entries.items():
            p = os.path.join(root, key)
            paths.append(p)
            ArchiveFileDB.to_list_path_str(val, p, paths)
        return list(paths)

    @staticmethod
    def to_nested_dict(path_str):
        if isinstance(path_str, str):
            path_str = path_str.split('/')
        if len(path_str) == 1:
            return {path_str[0]: None}
        else:
            pdict = {}
            pdict[path_str[0]] = ArchiveFileDB.to_nested_dict(path_str[1:])
            return pdict

    @staticmethod
    def append_data(entry, val):
        for k, v in entry.items():
            if v is None:
                entry[k] = val
            else:
                entry[k] = ArchiveFileDB.append_data(v, val)
        return entry

    @staticmethod
    def merge_dict(dict0, dict1):
        for k, v in dict1.items():
            if k in dict0:
                ArchiveFileDB.merge_dict(dict0[k], v)
            else:
                dict0[k] = v
        return dict0

    @property
    def mode(self):
        if 'b' not in self._mode:
            self._mode += 'b'
        return self._mode

    @mode.setter
    def mode(self, m):
        if self.mode == m:
            return
        self._mode = m
        if self._fileobj and not isinstance(self._fileobj, BytesIO):
            self._fileobj.close()

    @property
    def fileobj(self):
        if self._fileobj is None:
            mode = self.mode
            self._fileobj = open(self._filename, mode)
        return self._fileobj

    def get_docs(self, offset):
        """
        Provides an entry in the database.
        Arguments:
            offset: str to indicate the offset for unpacking the
                msgpack file or a string corresponding to the id of the entry
        """
        unpacker = msgpack.Unpacker(self.fileobj, raw=False)
        if isinstance(offset, str):
            self.fileobj.seek(0)
            sep = '%s_%s' % (self._sep, offset)
            offset = self.fileobj.read().find(sep.encode())
            if offset < 0:
                return
            offset += len(sep)
        self.fileobj.seek(offset)
        res = unpacker.unpack()
        return res

    @property
    def ids(self):
        if self._ids is None:
            self._ids = self.get_docs('IDS')
        return self._ids

    def _query(self, path_str):
        if path_str not in self.ids:
            return
        path_index = self.ids[path_str]
        data = self.get_docs(path_index)
        entry = ArchiveFileDB.to_nested_dict(path_str)
        if isinstance(data, dict):
            data = ArchiveFileDB.append_data(entry, list(data.values())[0])
            return data
        elif isinstance(data, list):
            d = {}
            for p in data:
                d = ArchiveFileDB.merge_dict(d, self._query(p))
            return d

    def query(self, entries, dtype='dict'):
        """
        Queries the database given a schema.
        Arguments:
            entries: A dictionary with a schema similar to the database
                entries but containing null values corresponding to the
                desired quantity.
            dtype: format of the outfile can be file, string, dict or metainfo
        """
        reduced_entries = self._reduce_to_section(entries)
        path_strs = ArchiveFileDB.to_list_path_str(reduced_entries, root='', paths=[])
        data = {}
        for path_str in path_strs:
            d = self._query(path_str)
            if d:
                data.update(d)
        res = {}
        if data:
            jdata = JSONdata(data)
            res = jdata.get_data(entries)
        if dtype == 'file':
            res = StringIO(json.dumps(res, indent=4))
        elif dtype == 'string':
            res = json.dumps(res)
        else:
            pass
        return res
