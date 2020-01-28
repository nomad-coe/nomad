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
Module for converting archive data into the new metainfo format.

In module ``ArchiveMetainfo``, the data is provided either from raw
json data or as a filename of an existing msgpack database. The metainfo
can then queried by providing a schema.

.. code-block: python
    am = ArchiveMetainfo("db.msg")
    for calc in am.calcs:
        c.section_run.section_single_configuration_calculation[0]({'energy_total':None})
"""

import numpy as np
from io import BytesIO

from nomad.metainfo import MSection, Quantity, SubSection
from nomad.metainfo.metainfo import MObjectMeta
from nomad.archive_library.filedb import ArchiveFileDB


class ArchiveMetainfo:
    """
    Converts archive data in json format to the new nomad metainfo model
    Arguments:
        archive_data: the archive data in json format or msgdb filename
        archive_schema: dict with the desired quantities as keys and None as placeholder
        for the values which are queried from the data
    """
    def __init__(self, archive_data, archive_schema=None):
        self._archive_data = archive_data
        self._archive_schema = archive_schema
        self.metainfo = None
        self._metacls = None
        self._calcs = {}
        self._calc_ids = []
        self._archive_db = None
        self._base_metacls = None
        self._base_metainfo = None
        self._base_data = None
        self._prefix = 'calc'
        self._load_archive_db()
        self._init_calcs()

    def _load_archive_db(self):
        if isinstance(self._archive_data, str):
            self._archive_db = ArchiveFileDB(self._archive_data)
        else:
            db = ArchiveFileDB(BytesIO(), mode='wb')
            if isinstance(self._archive_data, dict):
                for calc_id, run in self._archive_data.items():
                    db.add_data({calc_id: run})
            elif isinstance(self._archive_data, list):
                for entry in self._archive_data:
                    db.add_data(entry)
            db.create_db()
            self._archive_db = db

    def _init_calcs(self):
        for i in range(len(self.calc_ids)):
            calc_id = self.calc_ids[i]
            if self._archive_schema is None:
                self._calcs[calc_id] = self.base_metainfo
            else:
                data = self._archive_db.query({calc_id: self._archive_schema})[calc_id]
                self._calcs[calc_id] = self.base_metacls.m_from_dict(data)
            self._calcs[calc_id].archive_db = self._archive_db

    def __getitem__(self, key):
        if isinstance(key, str):
            calc = self._calcs.get(key, None)
            if calc:
                calc.calc_id = key
            return calc
        elif isinstance(key, int):
            calc_id = self._calc_ids[key]
            calc = self._calcs[calc_id]
            calc.calc_id = calc_id
            return calc
        else:
            calc_ids = self._calc_ids[key]
            calcs = []
            for calc_id in calc_ids:
                calc = self._calcs[calc_id]
                calc.calc_id = calc_id
                calcs.append(calc)
            return calcs

    def __len__(self):
        return len(self._calcs)

    def __iter__(self):
        self._n = -1
        return self

    def __next__(self):
        self._n += 1
        if self._n >= len(self):
            raise StopIteration
        calc = list(self._calcs.values())[self._n]
        calc.calc_id = list(self._calcs.keys())[self._n]
        calc.archive_db = self._archive_db
        return calc

    @staticmethod
    def get_path_from_section(content):
        path = content.m_path()
        path = path.split('/')
        s = ''
        for p in path:
            try:
                p = int(p)
                s += '[%s]' % p
            except ValueError:
                s += '/%s' % p
        return s[1:]

    @staticmethod
    def get_data_from_db(content, qschema):
        db = content.m_root().archive_db
        calc_id = content.m_root().calc_id
        root = calc_id + ArchiveMetainfo.get_path_from_section(content)
        qs = ArchiveFileDB.append_data(ArchiveFileDB.to_nested_dict(root), qschema)
        data = db.query(qs)
        return data

    @property
    def calcs(self):
        """
        Calculations in metainfo form which can be actively queried by using the get
        functionality and providing a schema
        """
        if not self._calcs:
            self._init_calcs()
        for calc_id, calc in self._calcs.items():
            calc.calc_id = calc_id
            calc.archive_db = self._archive_db
            yield calc

    @property
    def calc_ids(self):
        if not self._calc_ids:
            ids = self._archive_db.ids
            ids = [idx for idx in ids if '/' not in idx and idx]
            self._calc_ids = ids
        return self._calc_ids

    def _nullify_metainfo(self, metainfo):
        if hasattr(metainfo, 'm_contents'):
            for content in metainfo.m_contents():
                self._nullify_metainfo(content)
        return metainfo

    def _nullify_data(self, data):
        if not data:
            return
        elif isinstance(data, dict):
            for key, val in data.items():
                data[key] = self._nullify_data(val)
        elif isinstance(data, list) and isinstance(data[0], dict):
            for i in range(len(data)):
                data[i] = self._nullify_data(data[i])
        else:
            data = None
        return data

    @property
    def base_data(self):
        if self._base_data is None:
            calc_id = self.calc_ids[0]
            self._base_data = self._archive_db.query({calc_id: self._archive_schema})[calc_id]
        return self._base_data

    @property
    def base_metacls(self):
        """
        The base metaclass to apply a calculation
        """
        if self._base_metacls is None:
            name = self._prefix
            self._base_metacls = self._build_meta_cls(self.base_data, name)
        return self._base_metacls

    @property
    def base_metainfo(self):
        """
        The base metainfo to enable auto completion for each calc
        """
        if self._base_metainfo is None:
            metacls = self.base_metacls
            base_data = self._nullify_data(self.base_data)
            self._base_metainfo = metacls.m_from_dict(base_data)
        return self._base_metainfo

    def get_dtype(self, data):
        if isinstance(data, np.ndarray):
            if len(data) == 0:
                dtype = int
            else:
                dtype = self.get_dtype(data[0])
        else:
            dtype = type(data)
        return dtype

    def _to_meta_obj(self, content):
        if isinstance(content, Quantity):
            return content
        if isinstance(content, MObjectMeta):
            return SubSection(sub_section=content, repeats=content.repeats)
        else:
            if isinstance(content, list):
                content = np.array(content)
            dtype = self.get_dtype(content)
            if isinstance(content, np.ndarray):
                dtype = np.dtype(dtype)
                shape = np.shape(content)
                return Quantity(type=dtype, shape=shape)
            else:
                return Quantity(type=dtype)

    def _create_section(self, name, contents):
        contents['get'] = ArchiveMetainfo.get_data_from_db
        section = type(name.title(), (MSection,), contents)
        section.__call__ = ArchiveMetainfo.get_data_from_db
        return section

    def _build_meta_cls(self, data=None, name=None, return_section=True):
        if name is None:
            data = self._archive_data
            name = self._prefix
        if data is None:
            return
        if isinstance(data, dict):
            contents = {}
            for key, val in data.items():
                content = self._build_meta_cls(val, key, True)
                content = self._to_meta_obj(content)
                contents[key] = content
            if return_section:
                section = self._create_section(name, contents)
                section.repeats = False
                return section
            else:
                return contents

        elif isinstance(data, list):
            if not data:
                return self._to_meta_obj(data)
            if not isinstance(data[0], dict):
                return self._to_meta_obj(data)
            contents = {}
            for i in range(len(data)):
                content = self._build_meta_cls(data[i], name, False)
                contents.update(content)
            section = self._create_section(name, contents)
            section.repeats = True
            return section

        else:
            return self._to_meta_obj(data)

    def to_metainfo(self, data=None):
        if data is None:
            data = self._archive_data
        self.metainfo = self.base_metacls.m_from_dict(data)
