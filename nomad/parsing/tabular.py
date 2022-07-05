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

from typing import Union, List, Iterable, Dict, Callable, Set, Any, Tuple, cast
from memoization import cached
import os.path
import re


from nomad import utils
from nomad.units import ureg
from nomad.datamodel import EntryArchive, EntryData, ArchiveSection
from nomad.datamodel.context import Context
from nomad.metainfo import Section, Quantity, Package, Reference, SectionProxy, MSection, Property
from nomad.metainfo.metainfo import MetainfoError, SubSection
from nomad.parsing.parser import MatchingParser


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


class TableData(ArchiveSection):
    def normalize(self, archive, logger):
        super(TableData, self).normalize(archive, logger)

        for quantity in self.m_def.all_quantities.values():
            tabular_parser_annotation = quantity.m_annotations.get('tabular_parser', None)
            if tabular_parser_annotation:
                self.tabular_parser(quantity, archive, logger, **tabular_parser_annotation)

    def tabular_parser(self, quantity_def: Quantity, archive, logger, **kwargs):
        if not quantity_def.is_scalar:
            raise NotImplementedError('CSV parser is only implemented for single files.')

        value = self.m_get(quantity_def)
        if not value:
            return

        with archive.m_context.raw_file(self.data_file) as f:
            data = read_table_data(self.data_file, f, **kwargs)

        parse_columns(data, self)


class XLSOnlyTableData(ArchiveSection):
    def normalize(self, archive, logger):
        super(XLSOnlyTableData, self).normalize(archive, logger)

        for quantity in self.m_def.all_quantities.values():
            tabular_parser_annotation = quantity.m_annotations.get('tabular_parser', None)
            if tabular_parser_annotation:
                self.tabular_parser(quantity, archive, logger)

    def tabular_parser(self, quantity_def: Quantity, archive, logger):
        if not quantity_def.is_scalar:
            raise NotImplementedError('CSV parser is only implemented for single files.')

        value = self.m_get(quantity_def)
        if not value:
            return

        with archive.m_context.raw_file(self.data_file) as f:
            exlFile = XLSOnly_read_table_data(self.data_file, f)

        XLSOnly_parse_columns(exlFile, self)


m_package.__init_metainfo__()


@cached(max_size=10)
def _create_column_to_quantity_mapping(section_def: Section):
    mapping: Dict[str, Callable[[MSection, Any], MSection]] = {}

    def add_section_def(section_def: Section, path: List[Tuple[SubSection, Section]]):
        properties: Set[Property] = set()

        for quantity in section_def.all_quantities.values():
            if quantity in properties:
                continue
            properties.add(quantity)

            tabular_annotation = quantity.m_annotations.get('tabular', None)
            if tabular_annotation and 'name' in tabular_annotation:
                col_name = tabular_annotation['name']
            else:
                col_name = quantity.name
                if len(path) > 0:
                    col_name = f'{".".join([item[0].name for item in path])}.{col_name}'

            if col_name in mapping:
                raise MetainfoError(
                    f'The schema has non unique column names. {col_name} exists twice. '
                    f'Column names must be unique, to be used for tabular parsing.')

            def set_value(section: MSection, value, path=path, quantity=quantity, tabular_annotation=tabular_annotation):
                import numpy as np
                for sub_section, section_def in path:
                    next_section = section.m_get_sub_section(sub_section, -1)
                    if not next_section:
                        next_section = section_def.section_cls()
                        section.m_add_sub_section(sub_section, next_section, -1)
                    section = next_section

                if tabular_annotation and 'unit' in tabular_annotation:
                    value *= ureg(tabular_annotation['unit'])

                if isinstance(value, (int, float, str)):
                    value = np.array(value)

                if len(value.shape) == 1 and len(quantity.shape) == 0:
                    if len(value) == 1:
                        value = value[0]
                    elif len(value) == 0:
                        value = None
                    else:
                        raise MetainfoError(
                            'The shape of {quantity.name} does not match the given data.')
                elif len(value.shape) != len(quantity.shape):
                    raise MetainfoError(
                        'The shape of {quantity.name} does not match the given data.')

                section.m_set(quantity, value)

            mapping[col_name] = set_value

        for sub_section in section_def.all_sub_sections.values():
            if sub_section in properties or sub_section.repeats:
                continue
            next_base_section = sub_section.sub_section
            properties.add(sub_section)
            for sub_section_section in next_base_section.all_inheriting_sections + [next_base_section]:
                add_section_def(sub_section_section, path + [(sub_section, sub_section_section,)])

    add_section_def(section_def, [])
    return mapping


def parse_columns(pd_dataframe, section: MSection):
    '''
    Parses the given pandas dataframe and adds columns (all values as array) to
    the given section.
    '''
    import pandas as pd
    data: pd.DataFrame = pd_dataframe

    mapping = _create_column_to_quantity_mapping(section.m_def)  # type: ignore
    for column in data:
        if column in mapping:
            mapping[column](section, data.loc[:, column])


def XLSOnly_parse_columns(pd_exlFile, section: MSection):
    '''
    Parses the given pandas dataframe and adds columns (all values as array) to
    the given section.
    '''
    import pandas as pd
    exlFile: pd.ExcelFile = pd_exlFile

    mapping = _create_column_to_quantity_mapping(section.m_def)  # type: ignore
    for column in mapping:
        if '/' in column:
            sheet_name, col_name = column.split('/')
            data = pd.read_excel(exlFile, sheet_name=sheet_name, comment='#')
            if col_name in data:
                mapping[column](section, data.loc[:, col_name])
        else:
            data = pd.read_excel(exlFile, sheet_name=0, comment='#')
            if column in data:
                mapping[column](section, data.loc[:, column])


def parse_table(pd_dataframe, section_def: Section, logger):
    '''
    Parses the given pandas dataframe and creates a section based on the given
    section_def for each row. The sections are filled with the cells from
    their respective row.
    '''
    import pandas as pd
    data: pd.DataFrame = pd_dataframe
    sections: List[MSection] = []

    mapping = _create_column_to_quantity_mapping(section_def)  # type: ignore
    for row_index, row in data.iterrows():
        section = section_def.section_cls()
        try:
            for column in data:
                if column in mapping:
                    try:
                        mapping[column](section, row[column])
                    except Exception as e:
                        logger.error(
                            f'could not parse cell',
                            details=dict(row=row_index, column=column), exc_info=e)
        except Exception as e:
            logger.error(f'could not parse row', details=dict(row=row_index), exc_info=e)
        sections.append(section)

    return sections


def read_table_data(path, file_or_path=None, **kwargs):
    import pandas as pd
    if file_or_path is None:
        file_or_path = path
    if path.endswith('.xls') or path.endswith('.xlsx'):
        return pd.read_excel(
            file_or_path if isinstance(file_or_path, str) else file_or_path.name,
            **kwargs
        )
    else:
        return pd.read_csv(file_or_path, engine='python', **kwargs)


def XLSOnly_read_table_data(path, file_or_path=None, **kwargs):
    import pandas as pd
    if file_or_path is None:
        file_or_path = path
    if path.endswith('.xls') or path.endswith('.xlsx'):
        return pd.ExcelFile(file_or_path if isinstance(file_or_path, str) else file_or_path.name)


class TabularDataParser(MatchingParser):
    creates_children = True

    def __init__(self) -> None:
        super().__init__(
            name='parser/tabular', code_name='tabular data',
            mainfile_mime_re=r'text/.*|application/.*',
            mainfile_name_re=r'.*\.archive\.(csv|xlsx?)$')

    def _get_schema(self, filename: str, mainfile: str):
        dir = os.path.dirname(filename)
        match = re.match(r'^(.+\.)?([\w\-]+)\.archive\.(csv|xlsx?)$', os.path.basename(filename))
        if not match:
            return None

        schema_name = match.group(2)
        for extension in ['yaml', 'yml', 'json']:
            schema_file_base = f'{schema_name}.archive.{extension}'
            schema_file = os.path.join(dir, schema_file_base)
            if not os.path.exists(schema_file):
                continue
            return os.path.join(os.path.dirname(mainfile), schema_file_base)

        return None

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
            data = read_table_data(filename)
        except Exception:
            # If this cannot be parsed as a .csv file, we don't match with this file
            return False

        return [str(item) for item in range(0, data.shape[0])]

    def parse(
        self, mainfile: str, archive: EntryArchive, logger=None,
        child_archives: Dict[str, EntryArchive] = None
    ):
        if logger is None:
            logger = utils.get_logger(__name__)

        # We use mainfile to check the files existence in the overall fs,
        # and archive.metadata.mainfile to get an upload/raw relative schema_file
        schema_file = self._get_schema(mainfile, archive.metadata.mainfile)
        if schema_file is None:
            logger.error('Tabular data file without schema.', details=(
                'For a tabular file like name.schema.archive.csv, there has to be an '
                'uploaded schema like schema.archive.yaml'))
            return

        try:
            schema_archive = cast(Context, archive.m_context).load_raw_file(
                schema_file, archive.metadata.upload_id, None)
            package = schema_archive.definitions
            section_def = package.section_definitions[0]
        except Exception as e:
            logger.error('Could not load schema', exc_info=e)
            return

        if TableRow.m_def not in section_def.base_sections:
            logger.error('Schema for tabular data must inherit from TableRow.')
            return

        tabular_parser_annotation = section_def.m_annotations.get('tabular-parser', None)
        if tabular_parser_annotation:
            data = read_table_data(mainfile, **tabular_parser_annotation)
        else:
            data = read_table_data(mainfile)

        child_sections = parse_table(data, section_def, logger=logger)
        assert len(child_archives) == len(child_sections)

        table = Table()
        archive.data = table

        child_section_refs: List[MSection] = []
        for child_archive, child_section in zip(child_archives.values(), child_sections):
            child_archive.data = child_section
            child_section_refs.append(child_section)
            child_section.table_ref = table
        table.row_refs = child_section_refs
