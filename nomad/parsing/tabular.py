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

from typing import Union, List, Iterable, Dict
import csv

from nomad.datamodel.datamodel import EntryArchive, EntryData
from nomad.metainfo import Section, Quantity, Package, Reference, SectionProxy

from .parser import MatchingParser

# We define a simple base schema for tabular data. The parser will then generate more
# specialized sections based on the table headers. These specialized defintions will use
# this as base section definition.
# TODO maybe this should be moved to nomad.datamodel.metainfo
m_package = Package()


class TableRow(EntryData):
    ''' Represents the data in one row of a table. '''
    table_ref = Quantity(
        type=Reference(SectionProxy('Table')),
        description='A reference to the table that this row is contained in.')


class Table(EntryData):
    ''' Represents a table with many rows. '''
    row_refs = Quantity(
        type=Reference(TableRow.m_def), shape=['*'],
        description='References that connect to each row. Each row is stored in it individual entry.')


m_package.__init_metainfo__()


class TabularDataParser(MatchingParser):
    # TODO this is super simple and needs extension. Currently parses files like:
    #    header_0,header_1
    #    0_0,0_1
    #    1_0,1_1
    # TODO also extend tests/parsing/test_tabular
    def __init__(self) -> None:
        super().__init__(
            name='parser/tabular', code_name='tabular data',
            mainfile_name_re=r'.*\.csv$')

    def _read_cvs(self, filename: str) -> List[List[str]]:
        with open(filename, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='#')
            return [row for row in reader]

    def is_mainfile(
        self, filename: str, mime: str, buffer: bytes, decoded_buffer: str,
        compression: str = None
    ) -> Union[bool, Iterable[str]]:
        # We use the main file regex capabilities of the superclass to check if this is a
        # .csv file
        is_tabular = super().is_mainfile(filename, mime, buffer, decoded_buffer, compression)
        if not is_tabular:
            return False

        try:
            rows = self._read_cvs(filename)
        except Exception:
            # If this cannot be parsed as a .csv file, we don't match with this file
            return False

        return [str(item) for item in range(1, len(rows))]

    def parse(
        self, mainfile: str, archive: EntryArchive, logger=None,
        child_archives: Dict[str, EntryArchive] = None
    ):
        rows = self._read_cvs(mainfile)
        header = rows[0]
        rows = rows[1:]

        # Create a schema from the header and store it in the main archive.
        # We are creating a specialized TableRow section that has a quantity
        # for each column.
        # TODO do something smart with the headers to create more complex schemas.
        # e.g. parse headers like 'sec1/sec2/quantity_name?unit=m&type=float'
        table_row_def = Section(name='TableRow', base_sections=[TableRow.m_def])
        table_row_def.quantities = [
            Quantity(name=name, type=str) for name in header]
        archive.definitions = Package(section_definitions=[table_row_def])
        MyTableRow = table_row_def.section_cls

        # Create a Table as the main archive contents. This can hold references to the
        # rows.
        table = Table()
        archive.data = table
        table_rows = []

        # Create the data for each row by instantiating the generated schema and store it
        # in the child archives. Use references to connect with the Table section in the
        # main archive.
        for index, row in enumerate(rows):
            key = str(index + 1)
            archive = child_archives[key]
            archive.data = MyTableRow(table_ref=table, **dict(zip(header, row)))  # type: ignore
            table_rows.append(archive.data)
        table.row_refs = table_rows
