import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
)
from nomad.metainfo.legacy import LegacyDefinition

from nomad.datamodel.metainfo import general_experimental

m_package = Package(
    name='general_experimental_data_nomadmetainfo_json',
    description='None',
    a_legacy=LegacyDefinition(name='general.experimental.data.nomadmetainfo.json'))


class section_data(general_experimental.section_data):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_data'))

    entry_repository_url = Quantity(
        type=str,
        shape=[],
        description='''
        An URL to the entry on the repository, where the data is stored.
        ''',
        a_legacy=LegacyDefinition(name='entry_repository_url'))


m_package.__init_metainfo__()
