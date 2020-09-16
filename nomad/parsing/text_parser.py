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


import numpy as np
import re
import mmap
import logging


class Quantity:
    def __init__(self, name, re_pattern, str_operation=None, unit=None, dtype=None, comment=None):
        self.name = name
        self.re_pattern = re.compile(re_pattern.encode())
        self.unit = unit
        self.shape = None
        self.dtype = dtype
        self._value = None
        self._str_operation = str_operation
        self._comment = comment

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, val_in):

        def convert_value(val):
            if self._comment is not None:
                if val.strip()[0] == self._comment:
                    return

            if self._str_operation is not None:
                val = self._str_operation(val)

            else:
                val = val.strip().split() if isinstance(val, str) else val
                val = val[0] if len(val) == 1 else val

            def _convert(val):
                if isinstance(val, str):
                    if self.dtype is None:
                        if val.isdecimal():
                            val = int(val)
                        else:
                            try:
                                val = float(val)
                            except Exception:
                                pass

                    self.shape = []
                    return val

                elif type(val) in [np.ndarray, list]:
                    try:
                        dtype = float if self.dtype is None else self.dtype
                        val = np.array(val, dtype=dtype)
                        self.dtype = dtype
                        if np.all(np.mod(val, 1) == 0):
                            val = np.array(val, dtype=int)
                            self.dtype = int
                        self.shape = list(np.shape(val))

                    except Exception:
                        val = [_convert(v) for v in val]
                        self.dtype = None

                    return val

                else:
                    self.dtype = type(val)
                    self.shape = []
                    return val

            val = _convert(val)

            return val

        self._value = None if val_in is None else [convert_value(val) for val in val_in]

    def to_SI(self):
        pass


class UnstructuredTextFileParser:
    def __init__(self, mainfile, quantities, logger=None):
        self._mainfile = mainfile
        self.quantities = quantities
        self._file_mmap = None
        self.logger = logger if logger else logging

    @property
    def quantities(self):
        return self._quantities

    @quantities.setter
    def quantities(self, val):
        self._quantities = val
        self.quantities_mapping = {val[idx].name: idx for idx in range(len(val))}

    @property
    def file_mmap(self):
        if self._file_mmap is None:
            with open(self.mainfile) as f:
                self._file_mmap = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        return self._file_mmap

    @property
    def mainfile(self):
        return self._mainfile

    @mainfile.setter
    def mainfile(self, val):
        for quantity in self.quantities:
            quantity.value = None

        self._file_mmap = None
        self._mainfile = val

    def __getitem__(self, key):
        idx = self.quantities_mapping.get(key, None)
        if idx is None:
            return

        value = self.quantities[idx].value
        if value is None:
            self.parse(key)

        return self.quantities[idx].value

    def __setitem__(self, key, val):
        idx = self.quantities_mapping.get(key, None)
        if idx is None:
            return

        self.quantities[idx].value = val

    def keys(self):
        return self.quantities_mapping.keys()

    def items(self):
        for key in self.keys():
            yield key, self[key]

    def parse(self, key=None):
        if isinstance(key, str):
            key = [key]

        for quantity in self.quantities:
            if key is None or quantity.name in key:
                value = []
                for res in quantity.re_pattern.finditer(self.file_mmap):
                    value.append(''.join([group.decode() for group in res.groups() if group]))
                if not value:
                    continue

                try:
                    quantity.value = value
                except Exception:
                    self.logger.warn('Error setting value for %s ' % quantity.name)
                    pass
