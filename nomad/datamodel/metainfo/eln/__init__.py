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

from nomad import utils
from nomad.datamodel.data import EntryData
from nomad.datamodel.results import ELN, Results, Material
from nomad.metainfo import MSection, Package, Quantity, Datetime

m_package = Package(name='material_library')


class ElnBaseSection(MSection):
    '''
    A generic abstract base section for ELNs that provides a few commonly used properties.

    If you inherit from this section, but do not need some quantities, list those
    quantities in the `eln.hide` annotation of your inheriting section definition.

    Besides predefining some quantities, these base sections will add some metadata
    to NOMAD's search. A particular example are `tags`, if you define a string
    or enum quantity in your sections named `tags`, its values will be searchable.
    '''
    name = Quantity(
        type=str,
        description='A short human readable and descriptive name.',
        a_eln=dict(component='StringEditQuantity'))

    lab_id = Quantity(
        type=str,
        description='A id string that is unique at least for the lab that produced this data.',
        a_eln=dict(component='StringEditQuantity'))

    description = Quantity(
        type=str,
        description=(
            'A humand description. This provides room for human readable information '
            'that could not be captured in the ELN.'),
        a_eln=dict(component='RichTextEditQuantity'))

    def normalize(self, archive, logger):
        if isinstance(self, EntryData):
            if archive.data == self and self.name:
                archive.metadata.entry_name = self.name
            EntryData.normalize(self, archive, logger)

        if not archive.results:
            archive.results = Results(eln=ELN())
        if not archive.results.eln:
            archive.results.eln = ELN()

        for quantity in self.m_def.all_quantities.values():
            tabular_parser_annotation = quantity.m_annotations.get('tabular_parser', None)
            if tabular_parser_annotation:
                self.tabular_parser(quantity, archive, logger, **tabular_parser_annotation)

        if self.lab_id:
            if archive.results.eln.lab_ids is None:
                archive.results.eln.lab_ids = []
            archive.results.eln.lab_ids.append(self.lab_id)

        if getattr(self, 'name'):
            if archive.results.eln.names is None:
                archive.results.eln.names = []
            archive.results.eln.names.append(self.name)

        if getattr(self, 'description'):
            if archive.results.eln.descriptions is None:
                archive.results.eln.descriptions = []
            archive.results.eln.descriptions.append(self.description)

        if getattr(self, 'tags', None):
            if archive.results.eln.tags is None:
                archive.results.eln.tags = []
            tags = self.tags
            if isinstance(tags, list):
                archive.results.eln.tags.extend(tags)
            else:
                archive.results.eln.tags.append(tags)

        if not archive.results.eln.sections:
            archive.results.eln.sections = []
        archive.results.eln.sections.append(self.m_def.name)

    def tabular_parser(self, quantity_def: Quantity, archive, logger, **kwargs):
        from nomad.parsing.tabular import parse_columns, read_table_data

        if not quantity_def.is_scalar:
            raise NotImplementedError('CSV parser is only implemented for single files.')

        value = self.m_get(quantity_def)
        if not value:
            return

        with archive.m_context.raw_file(self.data_file) as f:
            data = read_table_data(self.data_file, f, **kwargs)

        parse_columns(data, self)


class ElnActivityBaseSecton(ElnBaseSection):
    '''
    A generic abstract base section for ELNs that provides a few commonly used for
    laboratory activities, e.g. processes, characterizations, measurements, etc.
    '''
    datetime = Quantity(
        type=Datetime,
        description='The date and time when this activity was done.',
        a_eln=dict(component='DateTimeEditQuantity'))

    method = Quantity(
        type=str,
        description='A short consistent handle for the applied method.')

    def normalize(self, archive, logger):
        super().normalize(archive, logger)

        if archive.results.eln.methods is None:
            archive.results.eln.methods = []
        if self.method:
            archive.results.eln.methods.append(self.method)
        else:
            archive.results.eln.methods.append(self.m_def.name)


class ElnWithFormulaBaseSection(ElnBaseSection):
    '''
    A generic abstract base section for ELNs that provides a few commonly used for
    items with a chemical formula, e.g. chemicals or samples.
    '''
    chemical_formula = Quantity(
        type=str,
        description=(
            'The chemical formula. This will be used directly and '
            'indirectly in the search. The formula will be used itself as well as '
            'the extracted chemical elements.'),
        a_eln=dict(component='StringEditQuantity'))

    def normalize(self, archive, logger):
        super().normalize(archive, logger)
        if logger is None:
            logger = utils.get_logger(__name__)
        from ase import Atoms
        if self.chemical_formula:
            if not archive.results:
                archive.results = Results()
            if not archive.results.material:
                archive.results.material = Material()
            material = archive.results.material

            try:
                atoms = Atoms(self.chemical_formula)
                material.elements = list(set(atoms.get_chemical_symbols()))
                material.chemical_formula_hill = atoms.get_chemical_formula(mode='hill')
                material.chemical_formula_reduced = atoms.get_chemical_formula(mode='reduce')
                material.chemical_formula_descriptive = self.chemical_formula
            except Exception as e:
                logger.warn('could not analyse chemical formula', exc_info=e)


class Chemical(ElnWithFormulaBaseSection):
    ''' A ELN base section that can be used for chemicals.'''
    pass


class Sample(ElnWithFormulaBaseSection):
    ''' A ELN base section that can be used for samples.'''
    pass


class Instrument(ElnBaseSection):
    ''' A ELN base section that can be used for instruments.'''
    def normalize(self, archive, logger):
        super().normalize(archive, logger)

        if self.name:
            if archive.results.eln.instruments is None:
                archive.results.eln.instruments = []
            archive.results.eln.instruments.append(self.name)


class Process(ElnActivityBaseSecton):
    ''' A ELN base section that can be used for processes.'''
    pass


class Measurement(ElnActivityBaseSecton):
    ''' A ELN base section that can be used for measurements.'''
    pass


m_package.__init_metainfo__()
