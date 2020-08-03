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

'''
The *parsing* module is an interface for the existing NOMAD-coe parsers.
This module redefines some of the old NOMAD-coe python-common functionality to create a
more coherent interface to the parsers.

Assumption about parsers
------------------------
For now, we make a few assumption about parsers
- they always work on the same *meta-info* version
- they have no conflicting python requirements
- they can be loaded at the same time and can be used within the same python process
- they are uniquely identified by a GIT URL and publicly accessible
- their version is uniquely identified by a GIT commit SHA

Each parser is defined via an instance of :class:`Parser`. The implementation :class:`LegacyParser` is used for most NOMAD-coe parsers.

.. autoclass:: nomad.parsing.Parser
    :members:

The are sub-classes for parsers with special purposes.

.. autoclass:: nomad.parsing.Parser
.. autoclass:: nomad.parsing.MatchingParser
.. autoclass:: nomad.parsing.MissingParser
.. autoclass:: nomad.parsing.BrokenParser
.. autoclass:: nomad.parsing.TemplateParser
.. autoclass:: nomad.parsing.GenerateRandomParser
.. autoclass:: nomad.parsing.ChaosParser
.. autoclass:: nomad.parsing.EmptyParser


The implementation :class:`LegacyParser` is used for most NOMAD-coe parsers.

.. autoclass:: nomad.parsing.LegacyParser


The parser definitions are available via the following two variables.

.. autodata:: nomad.parsing.parsers.parsers
.. autodata:: nomad.parsing.parsers.parser_dict

Parsers are reused for multiple calculations.

Parsers and calculation files are matched via regular expressions.

.. autofunction:: nomad.parsing.parsers.match_parser

Parsers in NOMAD-coe use a *backend* to create output. There are different NOMAD-coe
basends. In nomad@FAIRDI, we only currently only use a single backed. The following
classes provide a interface definition for *backends* as an ABC and a concrete implementation
based on nomad@fairdi's metainfo:

.. autoclass:: nomad.parsing.AbstractParserBackend
    :members:
.. autoclass:: nomad.parsing.Backend
    :members:
'''

from nomad.parsing.legacy import AbstractParserBackend, Backend, BackendError, LegacyParser
from nomad.parsing.parser import Parser, BrokenParser, MissingParser, MatchingParser
from nomad.parsing.artificial import TemplateParser, GenerateRandomParser, ChaosParser, EmptyParser
