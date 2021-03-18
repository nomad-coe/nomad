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
import re
import numpy as np
from xml.etree import ElementTree
from lxml import etree

from nomad.parsing.file_parser import FileParser


class XMLParser(FileParser):
    '''
    Parser for XML files using ElementTree.

    Arguments:
        mainfile: the file to be parsed
        logger: logger
        convert: specifies if quantities are converted automatically
    '''
    def __init__(self, mainfile: str = None, logger=None, **kwargs):
        super().__init__(mainfile, logger=logger, open=kwargs.get('open', None))
        self.convert = kwargs.get('convert', True)
        self.init_parameters()

    def init_parameters(self):
        '''
        Method to call after loading the xml file.
        '''
        self._elements = None

    @property
    def root(self):
        '''
        Returns the root of the XML tree.
        '''
        if self._file_handler is None:
            if self.mainfile is None:
                return
            try:
                self._file_handler = ElementTree.parse(self.mainfile_obj).getroot()
            except Exception:
                try:
                    # I cannot use the lxml XMLParser directly because it is not compatible with
                    # the ElementTree implementation.
                    xml = etree.parse(self.mainfile_mainfile_obj, parser=etree.XMLParser(recover=True))
                    self._file_handler = ElementTree.fromstring(etree.tostring(xml))
                except Exception:
                    pass

                self.logger.error('failed to load xml file')
                try:
                    # I cannot use the lxml XMLParser directly because it is not compatible with
                    # the ElementTree implementation.
                    xml = etree.parse(self.open(self.mainfile), parser=etree.XMLParser(recover=True))
                    self._file_handler = ElementTree.fromstring(etree.tostring(xml))
                except Exception:
                    pass
            self.init_parameters()

        return self._file_handler

    @property
    def elements(self):
        '''
        Returns a list of all elements in the XML tree.
        '''
        if self._elements is None:
            self._elements = self.root.findall('.//')

        return self._elements

    def parse(self, key, convert=None):
        '''
        Parse a quantity identified by key or an xpath-style path. Automatic conversion
        can be switch off by setting convert to False.
        '''
        _convert = convert if convert is not None else self._kwargs.get('convert', None)
        _convert = _convert if _convert is not None else self.convert
        if self._results is None:
            self._results = dict()

        if not self.root:
            return self

        key_in = key
        key = key.lstrip('/')
        if key.find('/') > 0:
            parent = os.path.dirname(key)
            child = os.path.basename(key)
            elements = self.root.findall(os.path.join('./', parent))
        else:
            elements = self.elements
            child = key

        val = []
        for element in elements:
            if child:
                v = element.attrib.get(child, None)
                if v is None:
                    v = element.findall(child)
                    v = [e.text if e.text is not None else '' for e in v]
                if v:
                    val.append(v)

            else:
                val.append(element.attrib)

        if not val:
            return self

        def convert_value(val_in):
            if isinstance(val_in, dict):
                val = dict()
                for k, v in val_in.items():
                    val[k] = convert_value(v)

            elif isinstance(val_in, str):
                # exponential formatting
                val = val_in.strip().split()
                if len(val) == 1:
                    val = val[0]
                    re_float = r'(\d*\.*\d*)d(\-*\+*\d+)'
                    val = re.sub(re_float, r'\1e\2', val)
                    if val.isdecimal():
                        val = int(val)
                    elif val == 'true' or val == 'false':
                        val = val == 'true'
                    else:
                        try:
                            val = float(val)
                        except Exception:
                            pass
                else:
                    val = [convert_value(v) for v in val]

            elif isinstance(val_in, list):
                try:
                    val = [v.split() if isinstance(v, str) else v for v in val_in]
                    val = [v[0] if (isinstance(v, list) and len(v) == 1) else v for v in val]
                    val = np.array(val, dtype=float)
                    if np.all(np.mod(val, 1) == 0):
                        val = np.array(val, dtype=int)
                except Exception:
                    val = [convert_value(v) for v in val_in]

            return val

        if _convert:
            val = convert_value(val)

        val = val[0] if len(val) == 1 else val

        self._results[key_in] = val
        return self
