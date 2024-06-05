import os.path
import re
from typing import Optional

import numpy as np
import yaml

NEXUS_AVAILABLE = True
try:
    from pynxtools.dataconverter import convert as pynxtools_converter
    from pynxtools.dataconverter import writer as pynxtools_writer
    from pynxtools.dataconverter.template import Template
    from pynxtools.definitions.dev_tools.utils.nxdl_utils import (
        get_app_defs_names,  # pylint: disable=import-error
    )
except ImportError:
    NEXUS_AVAILABLE = False
    pass

from nomad.datamodel.data import EntryData
from nomad.metainfo import MEnum, Package, Quantity
from nomad.units import ureg

m_package = Package(name='nexus_data_converter')


def create_eln_dict(archive):
    def transform(quantity_def, section, value, path):
        if quantity_def.unit:
            Q_ = ureg.Quantity
            val_unit = Q_(value, quantity_def.unit)

            default_display_unit = quantity_def.m_annotations.get(
                'eln', {'defaultDisplayUnit': None}
            ).defaultDisplayUnit
            if default_display_unit:
                val_unit = val_unit.to(default_display_unit)

            return dict(
                value=val_unit.magnitude.tolist()
                if isinstance(val_unit.magnitude, np.ndarray)
                else val_unit.magnitude,
                unit=str(format(val_unit.units, '~')),
            )
        return value

    def exclude(quantity_def, section):
        return quantity_def.name in ('reader', 'input_files', 'output', 'nxdl')

    eln_dict = archive.m_to_dict(transform=transform, exclude=exclude)
    del eln_dict['data']['m_def']

    return eln_dict


def write_yaml(archive, filename, eln_dict):
    with archive.m_context.raw_file(filename, 'w') as eln_file:
        yaml.dump(eln_dict['data'], eln_file, allow_unicode=True)


def populate_nexus_subsection(
    template: 'Template',
    app_def: str,
    archive,
    logger,
    output_file_path: Optional[str] = None,
    on_temp_file=False,
):
    """Populate nexus subsection in nomad from nexus template.

    There are three ways to populate nexus subsection from nexus template.
    1. First it writes a nexus file (.nxs), then the nexus subsectoin will be populated from
        that file.
    2. First it write the data in hdf5 datamodel (in a file in memory), later the nexus
        subsection will be populated from that in-memory file.
    3. (This is not yet done.) It directly poulate the nexus subsection from the template.

    Args:
        template: Nexus template.
        app_def: Name of application def NXxrd_pan.
        archive: AntryArchive section.
        output_file_path: Output file should be a relative path not absolute path.
        logger: nomad logger.
        on_temp_file: Whether data will be written in temporary disk, by default False.

    Raises:
        Exception: could not trigger processing from NexusParser
        Exception: could not trigger processing from NexusParser
    """
    _, nxdl_f_path = pynxtools_converter.helpers.get_nxdl_root_and_path(app_def)

    # Writing nxs file, parse and populate NeXus subsection:
    if output_file_path:
        archive.data.output = os.path.join(
            archive.m_context.raw_path(), output_file_path
        )
        pynxtools_writer.Writer(
            data=template, nxdl_f_path=nxdl_f_path, output_path=archive.data.output
        ).write()
        try:
            from nomad.parsing.nexus.nexus import NexusParser

            nexus_parser = NexusParser()
            nexus_parser.parse(
                mainfile=archive.data.output, archive=archive, logger=logger
            )
            try:
                archive.m_context.process_updated_raw_file(
                    output_file_path, allow_modify=True
                )
            except Exception as e:
                logger.error(
                    'could not trigger processing',
                    mainfile=archive.data.output,
                    exc_info=e,
                )
                raise e
            else:
                logger.info('triggered processing', mainfile=archive.data.output)
        except Exception as e:
            logger.error('could not trigger processing', exc_info=e)
            raise e

    # Write in temporary file and populate the NeXus section.
    elif not output_file_path or on_temp_file:
        output_file = 'temp_file.nxs'
        output_file = os.path.join(archive.m_context.raw_path(), output_file)
        logger.info(
            'No output NeXus file is found and data is being written temporary file.'
        )
        try:
            pynxtools_writer.Writer(
                data=template, nxdl_f_path=nxdl_f_path, output_path=output_file
            ).write()

            from nomad.parsing.nexus.nexus import NexusParser

            nexus_parser = NexusParser()
            nexus_parser.parse(mainfile=output_file, archive=archive, logger=logger)
            # Ensure no local reference with the hdf5file
        except Exception as e:
            logger.error('could not trigger processing', exc_info=e)
            raise e
        finally:
            if os.path.isfile(output_file):
                os.remove(output_file)


class ElnYamlConverter(EntryData):
    output = Quantity(
        type=str,
        description='Output yaml file to save all the data. Default: eln_data.yaml',
        a_eln=dict(component='StringEditQuantity'),
        a_browser=dict(adaptor='RawFileAdaptor'),
        default='eln_data.yaml',
    )

    def normalize(self, archive, logger):
        super(ElnYamlConverter, self).normalize(archive, logger)

        eln_dict = create_eln_dict(archive)
        write_yaml(archive, archive.data.output, eln_dict)


class NexusDataConverter(EntryData):
    reader = Quantity(
        type=MEnum(
            sorted(list(set(pynxtools_converter.get_names_of_all_readers())))
            if NEXUS_AVAILABLE
            else []
        ),
        description='The reader needed to run the Nexus converter.',
        a_eln=dict(component='AutocompleteEditQuantity'),
    )

    nxdl = Quantity(
        type=MEnum(sorted(list(set(get_app_defs_names()))) if NEXUS_AVAILABLE else []),
        description='The nxdl needed for running the Nexus converter.',
        a_eln=dict(component='AutocompleteEditQuantity'),
    )

    input_files = Quantity(
        type=str,
        shape=['*'],
        description='Input files needed to run the nexus converter.',
        a_eln=dict(component='FileEditQuantity'),
        a_browser=dict(adaptor='RawFileAdaptor'),
    )

    output = Quantity(
        type=str,
        description='Output Nexus filename to save all the data. Default: output.nxs',
        a_eln=dict(component='StringEditQuantity'),
        a_browser=dict(adaptor='RawFileAdaptor'),
        default='output.nxs',
    )

    def normalize(self, archive, logger):
        super(NexusDataConverter, self).normalize(archive, logger)

        raw_path = archive.m_context.raw_path()
        eln_dict = create_eln_dict(archive)

        if archive.data.input_files is None:
            archive.data.input_files = []

        if len(eln_dict['data']) > 0:
            write_yaml(archive, 'eln_data.yaml', eln_dict)

            if 'eln_data.yaml' not in archive.data.input_files:
                archive.data.input_files.append('eln_data.yaml')

        converter_params = {
            'reader': archive.data.reader,
            'nxdl': re.sub('.nxdl$', '', archive.data.nxdl),
            'input_file': [
                os.path.join(raw_path, file) for file in archive.data.input_files
            ],
            'output': os.path.join(raw_path, archive.data.output),
        }
        try:
            pynxtools_converter.logger = logger
            pynxtools_converter.helpers.logger = logger
            pynxtools_converter.convert(**converter_params)
        except Exception as e:
            logger.error(
                'could not convert to nxs', mainfile=archive.data.output, exc_info=e
            )
            raise e

        try:
            archive.m_context.process_updated_raw_file(
                archive.data.output, allow_modify=True
            )
        except Exception as e:
            logger.error(
                'could not trigger processing', mainfile=archive.data.output, exc_info=e
            )
            raise e
        else:
            logger.info('triggered processing', mainfile=archive.data.output)


m_package.__init_metainfo__()
