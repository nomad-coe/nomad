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
from nomad.metainfo import Package, Quantity, SubSection
from nomad.datamodel.data import EntryData
import nexusutils.dataconverter.readers.hall.reader as hall_reader  # pylint: disable=import-error
from .measurement import Measurement
from .hall_instrument import Instrument
from .nexus_to_msection import get_measurements, get_instrument


m_package = Package(name='ikz_hall')


class HallData(EntryData):
    """A parser for hall measurement data"""

    data_file = Quantity(
        type=str,
        a_eln=dict(component='FileEditQuantity'),
        a_browser=dict(adaptor='RawFileAdaptor'))

    measurements = SubSection(section_def=Measurement, repeats=True)

    def normalize(self, archive, logger):
        super(HallData, self).normalize(archive, logger)

        if not self.data_file:
            return

        logger.info('Parsing hall measurement measurement file.')
        with archive.m_context.raw_file(self.data_file, 'r', encoding='unicode_escape') as f:

            data_template = hall_reader.parse_txt(f.name)
            self.measurements = list(get_measurements(data_template))


class HallInstrument(EntryData):
    """Representation of an instrument"""

    data_file = Quantity(
        type=str,
        a_eln=dict(component='FileEditQuantity'),
        a_browser=dict(adaptor='RawFileAdaptor'))

    instrument = SubSection(section_def=Instrument)

    def normalize(self, archive, logger):
        super(HallInstrument, self).normalize(archive, logger)

        if not self.data_file:
            return

        logger.info('Parsing hall measurement instrument file.')
        with archive.m_context.raw_file(self.data_file, 'r', encoding='unicode_escape') as f:

            data_template = hall_reader.parse_txt(f.name)
            self.instrument = get_instrument(data_template, logger)


m_package.__init_metainfo__()
