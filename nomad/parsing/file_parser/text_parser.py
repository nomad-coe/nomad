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


import mmap
import io
import re
import numpy as np
import pint
from typing import List, Union, Callable, Type, Any

from nomad.parsing.file_parser import FileParser
from nomad.metainfo import Quantity as mQuantity
from nomad.utils import get_logger


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
                head,
                key,
                self._value,
                self._tail,
            )
        return self._re_pattern

    def __call__(self, text, repeats=True):
        values = []
        units = []
        if repeats:
            for res in self.re_pattern.finditer(text):
                unit = res.groupdict().get('__unit', None)
                values.append(
                    ''.join(
                        [
                            group.decode()
                            for group in res.groups()
                            if group and group != unit
                        ]
                    )
                )
                units.append(unit.decode() if unit is not None else None)
        else:
            res = self.re_pattern.search(text)
            if res is not None:
                unit = res.groupdict().get('__unit', None)
                units.append(unit.decode() if unit is not None else None)
                values.append(
                    ''.join(
                        [
                            group.decode()
                            for group in res.groups()
                            if group and group != unit
                        ]
                    )
                )


class Quantity:
    """
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

    """

    def __init__(
        self,
        quantity: Union[str, mQuantity],
        re_pattern: Union[str, ParsePattern],
        **kwargs,
    ):
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
            self.dtype = (
                quantity.type.type
                if isinstance(quantity.type, np.dtype)
                else quantity.type
            )
            self.unit = quantity.unit
            # check if metainfo shape has dependencies
            self.shape = quantity.shape
            if False in [str(i).isdigit() for i in self.shape]:
                self.shape = None
        # override metainfo
        self.dtype = kwargs.get('dtype', self.dtype)
        self.unit = kwargs.get('unit', self.unit)
        self.shape = kwargs.get('shape', self.shape)

        self._re_pattern: str = (
            re_pattern.re_pattern
            if isinstance(re_pattern, ParsePattern)
            else re_pattern
        )
        self.str_operation: Callable = kwargs.get('str_operation', None)
        self.sub_parser: TextParser = kwargs.get('sub_parser', None)
        self.repeats: bool = kwargs.get('repeats', False)
        self.convert: bool = kwargs.get('convert', True)
        self.flatten: bool = kwargs.get('flatten', True)
        self.reduce: bool = kwargs.get('reduce', True)
        self.comment: str = kwargs.get('comment', None)

    @property
    def re_pattern(self):
        """
        Returns a compiled re pattern.
        """
        if isinstance(self._re_pattern, str):
            re_pattern = self._re_pattern.replace('__unit', '__unit_%s' % self.name)
            self._re_pattern = re.compile(re_pattern.encode())
        return self._re_pattern

    @re_pattern.setter
    def re_pattern(self, val: str):
        self._re_pattern = val

    def to_data(self, val_raw: str):
        """
        Converts the parsed block into data.
        """

        def convert(val):
            if isinstance(val, str):
                if self.dtype is None:
                    if val.isdecimal():
                        return int(val)
                    else:
                        try:
                            return float(val)
                        except Exception:
                            pass
                else:
                    try:
                        return self.dtype(val)
                    except Exception:
                        pass

                return val

            elif isinstance(val, (list, np.ndarray)):
                try:
                    dtype = float if self.dtype is None else self.dtype
                    val_test = np.array(val, dtype=dtype)
                    if self.dtype is None:
                        if np.all(np.mod(val_test, 1) == 0):
                            val_test = np.array(val_test, dtype=int)
                            dtype = int
                    return val_test

                except Exception:
                    self.dtype = None
                    return [convert(v) for v in val]

            elif isinstance(val, dict):
                return {k: convert(v) for k, v in val.items()}

            else:
                return val

        if not val_raw:
            return

        if self.comment is not None:
            if val_raw.strip()[0] == self.comment:
                return

        data: Any = val_raw

        if self.str_operation is not None:
            data = self.str_operation(val_raw)

        elif self.flatten:
            data = val_raw.strip().split()
            if self.reduce:
                data = data[0] if len(data) == 1 else data

        if self.convert:
            data = convert(data)

        if isinstance(data, np.ndarray) and self.shape:
            try:
                data = np.reshape(data, self.shape)
            except Exception:
                pass

        return data

    def __repr__(self) -> str:
        if not self.sub_parser:
            return self.name
        sub_quantities = [q.name for q in self.sub_parser.quantities]
        return f'{self.name}({", ".join(sub_quantities[:5])}{"..." if len(sub_quantities) > 5 else ""})'


class TextParser(FileParser):
    """
    Parser for unstructured text files using the re module. The quantities to be parsed
    are given as a list of Quantity objects which specifies the regular expression. The mmap
    module is used to handle the file. By default, re.find_all is used to get matches
    for performance reasons. In this case, overlap is not tolerated in the re patterns.
    To avoid this, set findall to False to switch to re.finditer.

    Arguments:
        mainfile: the path to the file to be parsed
        quantities: list of Quantity objects to be parsed.
        logger: optional logger
        findall: if True will employ re.findall, otherwise re.finditer
        file_offset: offset in reading the file
        file_length: length of the chunk to be read from the file
    """

    def __init__(
        self,
        mainfile: str = None,
        quantities: List[Quantity] = None,
        logger=None,
        **kwargs,
    ):
        if logger is None:
            logger = get_logger(__name__)
        super().__init__(mainfile, logger=logger, open=kwargs.get('open', None))
        self._quantities: List[Quantity] = quantities
        self.findall: bool = kwargs.get('findall', True)
        self._kwargs = kwargs
        self._file_length: int = kwargs.get('file_length', 0)
        self._file_offset: int = kwargs.get('file_offset', 0)
        self._file_pad: int = 0
        if quantities is None:
            self.init_quantities()
        # check quantity patterns are valid
        re_has_group = re.compile(r'\(.+\)')
        for i in range(len(self._quantities) - 1, -1, -1):
            try:
                assert (
                    re_has_group.search(self._quantities[i].re_pattern.pattern.decode())
                    is not None
                )
            except Exception as e:
                self.logger.error(
                    'Invalid quantity pattern',
                    exc_info=e,
                    data=dict(quantity=self.quantities[i].name),
                )
                self._quantities.pop(i)
        self._re_findall: re.Pattern = None

    def copy(self):
        """
        Returns a copy of the object excluding the parsed results.
        """
        return TextParser(self.mainfile, self.quantities, self.logger, **self._kwargs)

    def init_quantities(self):
        """
        Initializes the quantities list.
        """
        self._quantities = []

    @property
    def quantities(self):
        """
        Returns the list of quantities to be parsed.
        """
        return self._quantities

    @quantities.setter
    def quantities(self, val: List[Quantity]):
        """
        Sets the quantities list.
        """
        self._file_handler = None
        self._results = None
        self._quantities = val

    @property
    def file_offset(self):
        """
        Integer offset in loading the file taking into account mmap pagination.
        """
        return self._file_offset

    @file_offset.setter
    def file_offset(self, val: int):
        """
        Sets starting point where the file is read.
        """
        self._file_pad = val % mmap.PAGESIZE
        self._file_offset = (val // mmap.PAGESIZE) * mmap.PAGESIZE
        self.reset()

    @property
    def file_length(self):
        """
        Length of the file chunk to be loaded.
        """
        return self._file_length

    @file_length.setter
    def file_length(self, val: int):
        """
        Sets the length of the file to be read.
        """
        self._file_length = val
        self.reset()

    @property
    def file_mmap(self):
        """
        Memory mapped representation of the file.
        """
        if self._file_handler is None:
            with self.open(self.mainfile) as f:
                if isinstance(f, io.TextIOWrapper):
                    self._file_handler = mmap.mmap(
                        f.fileno(),
                        self._file_length,
                        access=mmap.ACCESS_COPY,
                        offset=self._file_offset,
                    )
                    # set the extra chunk loaded before the intended offset to empty
                    self._file_handler[: self._file_pad] = b' ' * self._file_pad
                else:
                    self._file_handler = f.read()
            self._file_pad = 0
        return self._file_handler

    def keys(self):
        """
        Returns all the quantity names.
        """
        return [quantity.name for quantity in self.quantities]

    def items(self):
        """
        Returns an iterable name, value of the parsed quantities
        """
        for key in self.keys():
            yield key, self.get(key)

    def _add_value(self, quantity: Quantity, value: List[str], units):
        """
        Converts the list of parsed blocks into data and apply the corresponding units.
        """
        try:
            value_processed = [quantity.to_data(val) for val in value]
            for n, _ in enumerate(value_processed):
                unit = units[n] if units[n] else quantity.unit
                if not unit:
                    continue
                if isinstance(unit, str):
                    value_processed[n] = pint.Quantity(value_processed[n], unit)
                else:
                    value_processed[n] = value_processed[n] * unit

            if not quantity.repeats and value_processed:
                value_processed = value_processed[0]

            self._results[quantity.name] = value_processed
        except Exception:
            self.logger.warning(
                'Error setting value', data=dict(quantity=quantity.name)
            )

    def _parse_quantities(self, quantities: List[Quantity]):
        """
        Parse a list of quantities.
        """
        if len(self._results) == 0 and self._re_findall is not None:
            # attempt at optimization
            re_findall_b = self._re_findall
        else:
            re_findall = '|'.join([q.re_pattern.pattern.decode() for q in quantities])
            if len(quantities) == 1:
                # necessary to add a dummy variable to make multiple matches
                re_findall = '%s|(__dummy__)' % re_findall
            re_findall_b = re.compile(re_findall.encode())
            if self._re_findall is None:
                self._re_findall = re_findall_b

        # map matches to quantities
        matches = re.findall(re_findall_b, self.file_mmap)
        current_index = 0
        for quantity in quantities:
            values = []
            units = []
            n_groups = quantity.re_pattern.groups

            non_empty_matches = []
            for match in matches:
                non_empty_match = [
                    m for m in match[current_index : current_index + n_groups] if m
                ]
                if not non_empty_match:
                    continue
                non_empty_matches.append(non_empty_match)
            index_unit = quantity.re_pattern.groupindex.get(
                '__unit_%s' % quantity.name, None
            )
            for non_empty_match in non_empty_matches:
                try:
                    if index_unit is not None:
                        unit = non_empty_match.pop(index_unit - 1)
                        units.append(unit.decode())

                    else:
                        units.append(None)

                    values.append(' '.join([m.decode() for m in non_empty_match]))
                except Exception:
                    self.logger.error(
                        'Error parsing quantities.', data=dict(quantity=quantity.name)
                    )

            current_index += n_groups

            if not values:
                continue

            self._add_value(quantity, values, units)

    def _parse_quantity(self, quantity: Quantity):
        """
        Parse a single quantity.
        """
        value = []
        units = []
        re_matches = (
            quantity.re_pattern.finditer(self.file_mmap)
            if quantity.repeats
            else [quantity.re_pattern.search(self.file_mmap)]
        )
        for res in re_matches:
            if res is None:
                continue
            if quantity.sub_parser is not None:
                sub_parser = quantity.sub_parser.copy()
                sub_parser.mainfile = self.mainfile
                sub_parser.logger = self.logger
                sub_parser._file_handler = b' '.join([g for g in res.groups() if g])
                value.append(sub_parser.parse())

            else:
                try:
                    unit = res.groupdict().get('__unit_%s' % quantity.name, None)
                    units.append(unit.decode() if unit is not None else None)
                    value.append(
                        ' '.join(
                            [
                                group.decode()
                                for group in res.groups()
                                if group and group != unit
                            ]
                        )
                    )
                except Exception:
                    self.logger.error('Error parsing quantity.')

        if not value:
            return

        if quantity.sub_parser is not None:
            self._results[quantity.name] = value if quantity.repeats else value[0]

        else:
            self._add_value(quantity, value, units)

    def parse(self, key=None):
        """
        Triggers parsing of quantity with name key, if key is None will parse all quantities.

        Returns file parser.
        """
        if self._results is None:
            self._results = dict()

        if self.file_mmap is None:
            return self

        if self.findall:
            if len(self._results) > 1:
                return self

            n_results = 0
            while True:
                # use find all to parse quantities with no sub_parser.
                quantities_findall = [
                    q
                    for q in self.quantities
                    if q.name not in self._results and q.sub_parser is None
                ]
                if not quantities_findall:
                    break

                # recursively parse quantities
                self._parse_quantities(quantities_findall)

                if n_results == len(self._results):
                    # will stop if no more matches are found
                    break
                n_results = len(self._results)

            for quantity in self._quantities:
                if quantity.sub_parser is not None:
                    self._parse_quantity(quantity)

        else:
            for quantity in self._quantities:
                if quantity.name == key or key is None:
                    if quantity.name not in self._results:
                        self._parse_quantity(quantity)

        # free up memory
        if isinstance(self._file_handler, mmap.mmap) and self.findall:
            self._file_handler.close()
            self._file_handler = b' '

        return self

    def clear(self):
        """
        Deletes the file mapping for all sub parsers.
        """
        for quantity in self.quantities:
            if quantity.sub_parser is not None:
                quantity.sub_parser.clear()
        self._file_handler = None


class DataTextParser(TextParser):
    """
    Parser for structured data text files using numpy.loadtxt

    Arguments:
        mainfile: the file to be parsed
        dtype: data type
    """

    def __init__(self, **kwargs):
        self._dtype: Type = kwargs.get('dtype', float)
        self._mainfile_contents: str = kwargs.get('mainfile_contents', '')
        super().__init__(**kwargs)

    def parse(self, key=None):
        super().parse(key=key)
        if key == 'data':
            try:
                data = None
                if self.mainfile is not None:
                    data = np.loadtxt(self.mainfile)
                else:
                    if not self._mainfile_contents and self.mainfile_obj:
                        with self.open_mainfile_obj() as mainfile_obj:
                            self._mainfile_contents = mainfile_obj.read()
                    if self._mainfile_contents:
                        buffer = self._mainfile_contents
                        if isinstance(buffer, str):
                            buffer = buffer.encode()
                        if buffer:
                            data = np.frombuffer(buffer, dtype=self._dtype)
                if data is not None:
                    self._results['data'] = data
            except Exception:
                self.logger.error('Failed to load data file.')
