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


import logging
import mmap
import io
import re
import numpy as np
import pint
from typing import List, Union, Callable, Type

from nomad.parsing.file_parser import FileParser
from nomad.metainfo import Quantity as mQuantity


class ParsePattern:
    def __init__(self, **kwargs):
        self._head = kwargs.get('head', '')
        self._key = kwargs.get('key', '')
        value = kwargs.get('value', 're_float_array')
        if value.startswith('re_'):
            token = ''
            if 'float' in value:
                token += r'Ee\+\d\.\-'
            if 'int' in value:
                token += r'\d'
            if 'str' in value:
                token += r'\w'
            if 'array' in value:
                token += r' '
            value = r'[%s]+' % token
        self._value = value
        self._tail = kwargs.get('tail', '\n')
        self._re_pattern = None

    @property
    def re_pattern(self):
        if self._re_pattern is None:
            head = r'%s[\s\S]*?' % self._head if self._head else ''
            key = r'%s\s*\:*\=*\s*' % self._key if self._key else ''
            self._re_pattern = r'%s%s\s*\:*\=*\s*(%s)%s' % (
                head, key, self._value, self._tail)
        return self._re_pattern

    def __call__(self, text, repeats=True):
        values = []
        units = []
        if repeats:
            for res in self.re_pattern.finditer(text):
                unit = res.groupdict().get('__unit', None)
                values.append(
                    ''.join([group.decode() for group in res.groups() if group and group != unit]))
                units.append(unit.decode() if unit is not None else None)
        else:
            res = self.re_pattern.search(text)
            if res is not None:
                unit = res.groupdict().get('__unit', None)
                units.append(unit.decode() if unit is not None else None)
                values.append(''.join(
                    [group.decode() for group in res.groups() if group and group != unit]))


class Quantity:
    '''
    Class to define a quantity to be parsed in the TextParser.

    Arguments:
        quantity: string to identify the name or a metainfo quantity to initialize the
            quantity object.
        re_pattern: pattern to be used by re for matching. Ideally, overlaps among
            quantities for a given parser should be avoided.
        sub_parser: instance of TextParser to perform local parsing
            within a matched block
        str_operation: external function to be performed on a matched block
        dtype: data type of the quantity
        unit: unit of the quantity
        shape: shape of the quantity
        repeats: denotes if multiple matches are expected
        convert: switch automatic data type conversion
        comment: character to denote a line to be ignored

    '''
    def __init__(self, quantity: Union[str, mQuantity], re_pattern: Union[str, ParsePattern], **kwargs):
        self.name: str
        self.dtype: str
        self.unit: str
        self.shape: List[int]
        if isinstance(quantity, str):
            self.name = quantity
            self.dtype = None
            self.unit = None
            self.shape = None
        elif isinstance(quantity, mQuantity):
            self.name = quantity.name
            self.dtype = quantity.type.type if isinstance(quantity.type, np.dtype) else quantity.type
            self.unit = quantity.unit
            # check if metainfo shape has dependencies
            self.shape = quantity.shape
            if False in [str(i).isdigit() for i in self.shape]:
                self.shape = None
        # override metainfo
        self.dtype = kwargs.get('dtype', self.dtype)
        self.unit = kwargs.get('unit', self.unit)
        self.shape = kwargs.get('shape', self.shape)

        self._re_pattern: str = re_pattern.re_pattern if isinstance(
            re_pattern, ParsePattern) else re_pattern
        self._str_operation: Callable = kwargs.get('str_operation', None)
        self._sub_parser: TextParser = kwargs.get('sub_parser', None)
        self.repeats: bool = kwargs.get('repeats', False)
        self.convert: bool = kwargs.get('convert', True)
        self.flatten: bool = kwargs.get('flatten', True)
        self.comment: str = kwargs.get('comment', None)

    @property
    def re_pattern(self):
        '''
        Returns a compiled re pattern.
        '''
        if isinstance(self._re_pattern, str):
            re_pattern = self._re_pattern.replace('__unit', '__unit_%s' % self.name)
            self._re_pattern = re.compile(re_pattern.encode())
        return self._re_pattern

    @re_pattern.setter
    def re_pattern(self, val: str):
        self._re_pattern = val

    @property
    def str_operation(self):
        return self._str_operation

    @str_operation.setter
    def str_operation(self, val: Callable):
        self._str_operation = val

    def to_data(self, val_in: List[str]):
        '''
        Converts the parsed block into data.
        '''
        def process(val):
            if self.comment is not None:
                if val.strip()[0] == self.comment:
                    return

            if self.str_operation is not None:
                val = self.str_operation(val)

            elif self.flatten:
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
                    else:
                        try:
                            val = self.dtype(val)
                        except Exception:
                            pass

                    return val

                elif type(val) in [np.ndarray, list]:
                    try:
                        dtype = float if self.dtype is None else self.dtype
                        val_test = np.array(val, dtype=dtype)
                        if self.dtype is None:
                            if np.all(np.mod(val_test, 1) == 0):
                                val_test = np.array(val_test, dtype=int)
                        val = val_test

                    except Exception:
                        val = [_convert(v) for v in val]

                    return val

                elif isinstance(val, dict):
                    for k, v in val.items():
                        val[k] = _convert(v)
                    return val

                else:
                    return val

            if self.convert:
                val = _convert(val)

            if isinstance(val, np.ndarray) and self.shape:
                try:
                    val = np.reshape(val, self.shape)
                except Exception:
                    pass

            return val

        val_out = [process(val) for val in val_in]

        if isinstance(val_out[0], np.ndarray):
            self.dtype = val_out[0].dtype  # type: ignore

        return val_out


class DataTextParser(FileParser):
    '''
    Parser for structured data text files using numpy.loadtxt

    Arguments:
        mainfile: the file to be parsed
        dtype: data type
    '''
    def __init__(self, **kwargs):
        self._dtype: Type = kwargs.get('dtype', float)
        mainfile: str = kwargs.get('mainfile', None)
        self._mainfile_contents: str = kwargs.get('mainfile_contents', None)
        logger = kwargs.get('logger', None)
        logger = logger if logger is not None else logging
        super().__init__(mainfile, logger=logger, open=kwargs.get('open', None))
        self.init_parameters()

    def init_parameters(self):
        '''
        Method to call after loading data.
        '''
        pass

    @property
    def data(self):
        '''
        Returns the loaded data
        '''
        if self._file_handler is None:
            try:
                if self.mainfile is not None:
                    self._file_handler = np.loadtxt(self.mainfile)
                else:
                    if self._mainfile_contents is None:
                        self._mainfile_contents = self.mainfile_obj.read()
                    self._file_handler = np.frombuffer(self._mainfile_contents, dtype=self._dtype)
            except Exception:
                return

            self.init_parameters()
        return self._file_handler


class TextParser(FileParser):
    '''
    Parser for unstructured text files using the re module. The quantities to be parsed
    are given as a list of Quantity objects which specifies the re pattern. The mmap
    module is used to handle the file. By default, re.find_all is used to get matches
    for performance reasons. In this case, overlap is not tolerated in the re patterns.
    To avoid this, set findall to False to switch to re.finditer.

    Arguments:
        mainfile: the file to be parsed
        quantities: list of Quantity objects to be parsed.
        logger: optional logger
        findall: switches between using re.findall and re.finditer
        file_offset: offset in reading the file
        file_length: length of the chunk to be read from the file
    '''
    def __init__(self, mainfile=None, quantities=None, logger=None, findall=True, **kwargs):
        super().__init__(mainfile, logger=logger, open=kwargs.get('open', None))
        self._quantities: List[Quantity] = quantities
        self.findall: bool = findall
        self._kwargs = kwargs
        self._file_length: int = kwargs.get('file_length', 0)
        self._file_offset: int = kwargs.get('file_offset', 0)
        self._file_pad: int = 0
        if quantities is None:
            self.init_quantities()
        # check quantity patterns are valid
        re_has_group = re.compile(r'\(.+\)')
        for i in range(len(self.quantities) - 1, -1, -1):
            valid = True
            try:
                if re_has_group.search(self.quantities[i].re_pattern.pattern.decode()) is None:
                    valid = False
            except Exception:
                valid = False
            if not valid:
                self.logger.error(
                    'Invalid quantity pattern', data=dict(quantity=self.quantities[i].name))
                self.quantities.pop(i)
        self._re_findall = None

    def copy(self):
        '''
        Returns a copy of the object excluding the parsed results.
        '''
        return TextParser(
            self.mainfile, self.quantities, self.logger, **self._kwargs)

    def init_quantities(self):
        '''
        Initializes the quantities list.
        '''
        self._quantities = []

    @property
    def quantities(self):
        return self._quantities

    @quantities.setter
    def quantities(self, val):
        self._quantities = val

    @property
    def file_offset(self):
        '''
        Integer offset in loading the file taking into account mmap pagination.
        '''
        return self._file_offset

    @file_offset.setter
    def file_offset(self, val):
        self._file_pad = val % mmap.PAGESIZE
        self._file_offset = (val // mmap.PAGESIZE) * mmap.PAGESIZE

    @property
    def file_length(self):
        '''
        Length of the file chunk to be loaded.
        '''
        return self._file_length

    @file_length.setter
    def file_length(self, val):
        self._file_length = val

    @property
    def file_mmap(self):
        '''
        Memory mapped representation of the file.
        '''
        if self._file_handler is None:
            with self.open(self.mainfile) as f:
                if isinstance(f, io.TextIOWrapper):
                    self._file_handler = mmap.mmap(
                        f.fileno(), self._file_length, access=mmap.ACCESS_COPY,
                        offset=self._file_offset)
                    # set the extra chunk loaded before the intended offset to empty
                    self._file_handler[:self._file_pad] = b' ' * self._file_pad
                else:
                    self._file_handler = f.read()
            self._file_pad = 0
        return self._file_handler

    def keys(self):
        '''
        Returns all the quantity names.
        '''
        return [quantity.name for quantity in self.quantities]

    def items(self):
        '''
        Returns an iterable name, value of the parsed quantities
        '''
        for key in self.keys():
            yield key, self.get(key)

    def _parse_quantities(self, quantities):
        if len(self._results) == 0 and self._re_findall is not None:
            # maybe an opt
            re_findall = self._re_findall
        else:
            re_findall = '|'.join([q.re_pattern.pattern.decode() for q in quantities])
            if len(quantities) == 1:
                # necessary to add a dummy variable to make multiple matches
                re_findall = '%s|(__dummy__)' % re_findall
            re_findall = re_findall.encode()
            if self._re_findall is None:
                self._re_findall = re.compile(re_findall)

        # map matches to quantities
        matches = re.findall(re_findall, self.file_mmap)
        current_index = 0
        for i in range(len(quantities)):
            values = []
            units = []
            n_groups = quantities[i].re_pattern.groups

            non_empty_matches = []
            for match in matches:
                non_empty_match = [m for m in match[current_index: current_index + n_groups] if m]
                if not non_empty_match:
                    continue
                non_empty_matches.append(non_empty_match)
            index_unit = quantities[i].re_pattern.groupindex.get(
                '__unit_%s' % quantities[i].name, None)
            for non_empty_match in non_empty_matches:
                if index_unit is not None:
                    unit = non_empty_match.pop(index_unit - 1)
                    units.append(unit.decode())

                else:
                    units.append(None)

                values.append(' '.join([m.decode() for m in non_empty_match]))

            current_index += n_groups

            if not values:
                continue

            try:
                value_processed = quantities[i].to_data(values)
                for j in range(len(value_processed)):
                    unit = units[j] if units[j] else quantities[i].unit
                    if not unit:
                        continue
                    if isinstance(unit, str):
                        value_processed[j] = pint.Quantity(value_processed[j], unit)
                    else:
                        value_processed[j] = value_processed[j] * unit

                if not quantities[i].repeats and value_processed:
                    value_processed = value_processed[0]

                self._results[quantities[i].name] = value_processed

            except Exception:
                self.logger.warn('Error setting value', data=dict(quantity=quantities[i].name))
                pass

    def _parse_quantity(self, quantity):

        value = []
        units = []
        re_matches = quantity.re_pattern.finditer(self.file_mmap) if quantity.repeats else [
            quantity.re_pattern.search(self.file_mmap)]
        for res in re_matches:
            if res is None:
                continue
            if quantity._sub_parser is not None:
                span = np.array(res.span()) + self.file_offset
                sub_parser = quantity._sub_parser.copy()
                sub_parser.mainfile = self.mainfile
                sub_parser.logger = self.logger
                if (span[1] - span[0]) < mmap.PAGESIZE or True:
                    # self.logger.warn(
                    #     'Cannot use sub parser on quantity %s with blocks with size <'
                    #     '%d. Will try to parse string' % (quantity.name, mmap.PAGESIZE))
                    sub_parser._file_handler = b' '.join([g for g in res.groups() if g])
                else:
                    sub_parser.file_offset = span[0]
                    sub_parser.file_length = span[1] - sub_parser.file_offset
                value.append(sub_parser.parse())

            else:
                unit = res.groupdict().get('__unit_%s' % quantity.name, None)
                units.append(unit.decode() if unit is not None else None)
                value.append(' '.join(
                    [group.decode() for group in res.groups() if group and group != unit]))

        if not value:
            return

        if quantity._sub_parser is not None:
            self._results[quantity.name] = value if quantity.repeats else value[0]

        else:
            try:
                value_processed = quantity.to_data(value)
                for i in range(len(value_processed)):
                    unit = units[i] if units[i] else quantity.unit
                    if not unit:
                        continue
                    if isinstance(unit, str):
                        value_processed[i] = pint.Quantity(value_processed[i], unit)
                    else:
                        value_processed[i] = value_processed[i] * unit

                if not quantity.repeats and value_processed:
                    value_processed = value_processed[0]

                self._results[quantity.name] = value_processed
            except Exception:
                self.logger.warn('Error setting value', data=dict(quantity=quantity.name))
                pass

    def parse(self, key=None):
        '''
        Triggers parsing of all quantities if key is not provided.
        '''
        if self._results is None:
            self._results = dict()

        if self.file_mmap is None:
            return self

        if self.findall:
            if len(self._results) > 1:
                return self

            n_results = 0
            while True:
                quantities_findall = [
                    q for q in self.quantities if q.name not in self._results and q._sub_parser is None]
                if not quantities_findall:
                    break

                # recursively parse quantities
                self._parse_quantities(quantities_findall)

                if n_results == len(self._results):
                    break
                n_results = len(self._results)

            for quantity in self._quantities:
                if quantity._sub_parser is not None:
                    self._parse_quantity(quantity)

        else:
            for quantity in self._quantities:
                if quantity.name == key or key is None:
                    if quantity.name not in self._results:
                        self._parse_quantity(quantity)

        return self
