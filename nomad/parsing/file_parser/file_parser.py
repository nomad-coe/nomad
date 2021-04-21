# Copyright 2018 Markus Scheidgen
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


import os
import logging
import pint
from typing import Any, Dict
import gzip
import bz2
import lzma
import tarfile


class FileParser:
    '''
    Base class for parsers. The parse method specific to a file type
    should be implemented in the corresponding child class. The parsed quantities are
    stored in results. One can access a quantity by using the get method or as attribute.

    Arguments:
        mainfile: the file to be parsed
        logger: optional logger
    '''
    def __init__(self, mainfile=None, logger=None, open=None):
        self._mainfile: Any = None
        self._mainfile_obj: Any = None
        if isinstance(mainfile, str):
            self._mainfile = os.path.abspath(mainfile)
            self._mainfile_obj = None
        elif hasattr(mainfile, 'name'):
            self._mainfile = mainfile.name
            self._mainfile_obj = mainfile
        self._open = open
        self.logger = logger if logger is not None else logging
        self._results: Dict[str, Any] = None
        # a key is necessary for xml parsers, where parsing is done dynamically
        self._key: str = None
        self._kwargs: Dict[str, Any] = dict()
        self._file_handler: Any = None

    def init_parameters(self):
        pass

    @property
    def results(self):
        if self._results is None:
            self._results = dict()
        if self._key not in self._results:
            self.parse(self._key, **self._kwargs)

        return self._results

    @property
    def maindir(self):
        return os.path.dirname(self._mainfile)

    @property
    def mainfile_obj(self):
        if self._mainfile_obj is None:
            try:
                self._mainfile_obj = self.open(self._mainfile)
            except Exception:
                pass

        return self._mainfile_obj

    @property
    def mainfile(self):
        if self._mainfile is None:
            return

        if self._mainfile_obj is None and not os.path.isfile(self._mainfile):
            return
        return self._mainfile

    @mainfile.setter
    def mainfile(self, val):
        self._results = None
        self._file_handler = None
        self._mainfile = None
        if isinstance(val, str):
            self._mainfile = os.path.abspath(val)
            self._mainfile_obj = None
        elif hasattr(val, 'name'):
            self._mainfile = val.name
            self._mainfile_obj = val
        self.init_parameters()

    def open(self, mainfile):
        open_file = self._open
        if open_file is None:
            if mainfile.endswith('.gz'):
                open_file = gzip.open
            elif mainfile.endswith('.bz2'):
                open_file = bz2.open
            elif mainfile.endswith('.xz'):
                open_file = lzma.open
            elif mainfile.endswith('.tar'):
                open_file = tarfile.open
            else:
                open_file = open
        return open_file(mainfile)

    def get(self, key: str, default: Any = None, unit: Any = None, **kwargs):
        '''
        Returns the parsed result for quantity with name key. If quantity is not in
        results default will be returned. A pint unit can be provided which is attached
        to the returned value.
        '''
        if self.mainfile is None:
            return default

        self._key = key
        self._kwargs = kwargs
        val = self.results.get(key, None)
        if val is None:
            val = default

        if val is None:
            return

        if unit is not None:
            if isinstance(unit, pint.Quantity) or isinstance(unit, pint.Unit):
                val = val * unit

            elif isinstance(val, pint.Quantity):
                val = val.to(unit)

            else:
                val = pint.Quantity(val, unit)

        return val

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.get(key)
        elif isinstance(key, int):
            return self[int]

    def __getattr__(self, key):
        if self._results is None:
            self._results = dict()

        return self._results.get(key, None)

    def parse(self, quantity_key: str = None, **kwargs):
        pass
