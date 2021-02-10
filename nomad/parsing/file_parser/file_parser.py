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


class FileParser:
    '''
    Base class for parsers. The parse method implemented here simply sets the parsed
    quantities as attributes of the class. The parse method specific to a file type
    should be implemented in the corresponding child class. The parsed quantities are
    stored in results. One can access a quantity by using the get method.

    Arguments:
        mainfile: the file to be parsed
        logger: optional logger
    '''
    def __init__(self, mainfile: str, logger=None):
        self._mainfile: str = os.path.abspath(mainfile) if mainfile else mainfile
        self.logger = logger if logger else logging
        self._results: Dict[str, Any] = None
        # a key is necessary for xml parsers, where parsing is done dynamically
        self._key: str = None
        self._kwargs: Dict[str, Any] = None
        self._file_handler: Any = None

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
    def mainfile(self):
        if self._mainfile is None:
            return

        if not os.path.isfile(self._mainfile):
            return
        return self._mainfile

    @mainfile.setter
    def mainfile(self, val):
        self._results = None
        self._file_handler = None
        self._mainfile = os.path.abspath(val) if val is not None else val

    @property
    def open(self):
        if self.mainfile.endswith('.gz'):
            open_file = gzip.open
        elif self.mainfile.endswith('.bz2'):
            open_file = bz2.open
        elif self.mainfile.endswith('.xz'):
            open_file = lzma.open
        else:
            open_file = open
        return open_file

    def get(self, key: str, default: Any = None, unit: str = None, **kwargs):
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
            if isinstance(unit, pint.Quantity):
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

    def parse(self, quantity_key: str = None):
        '''
        Sets quantities in result as class attributes.
        '''
        for key, val in self._results.items():
            try:
                setattr(self, key, val)
            except Exception:
                pass
        return self
