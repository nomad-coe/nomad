import re
import yaml
import os.path

from nomad.datamodel.data import EntryData
from nomad.metainfo import Package, Quantity, MEnum

from nexusutils.dataconverter.convert import get_names_of_all_readers, convert
from nexusutils.nexus.nexus import get_app_defs_names


m_package = Package(name='nexus_data_converter')


class NexusDataConverter(EntryData):

    reader = Quantity(
        type=MEnum(get_names_of_all_readers()),
        description='The reader needed to run the Nexus converter.',
        a_eln=dict(component='AutocompleteEditQuantity'))

    nxdl = Quantity(
        type=MEnum(get_app_defs_names()),
        description='The nxdl needed for running the Nexus converter.',
        a_eln=dict(component='AutocompleteEditQuantity'))

    input_files = Quantity(
        type=str,
        shape=['*'],
        description='Input files needed to run the nexus converter.',
        a_eln=dict(component='FileEditQuantity'),
        a_browser=dict(adaptor='RawFileAdaptor'))

    output = Quantity(
        type=str,
        description='Output Nexus filename to save all the data. Default: output.nxs',
        a_eln=dict(component='StringEditQuantity'),
        a_browser=dict(adaptor='RawFileAdaptor'),
        default='output.nxs')

    def normalize(self, archive, logger):
        super(NexusDataConverter, self).normalize(archive, logger)

        raw_path = archive.m_context.raw_path()

        def transform(quantity_def, section, value, path):
            if quantity_def.unit:
                return dict(value=value, unit=str(format(quantity_def.unit, '~')))
            return value

        def exclude(quantity_def, section):
            return quantity_def.name in ('reader', 'input_files', 'output', 'nxdl')

        eln_dict = archive.m_to_dict(transform=transform, exclude=exclude)
        del eln_dict['data']['m_def']
        with archive.m_context.raw_file('eln_data.yaml', 'w') as eln_file:
            yaml.dump(eln_dict['data'], eln_file, allow_unicode=True)

        if archive.data.input_files is None:
            archive.data.input_files = []

        converter_params = {
            'reader': archive.data.reader,
            'nxdl': re.sub('.nxdl$', '', archive.data.nxdl),
            'input_file': [
                os.path.join(raw_path, file)
                for file in archive.data.input_files],
            'output': os.path.join(raw_path, archive.data.output)
        }
        try:
            convert(**converter_params)
        except Exception as e:
            logger.error('could not convert to nxs', mainfile=archive.data.output, exc_info=e)
            raise e

        try:
            archive.m_context.process_updated_raw_file(archive.data.output, allow_modify=True)
        except Exception as e:
            logger.error('could not trigger processing', mainfile=archive.data.output, exc_info=e)
            raise e
        else:
            logger.info('triggered processing', mainfile=archive.data.output)


m_package.__init_metainfo__()
