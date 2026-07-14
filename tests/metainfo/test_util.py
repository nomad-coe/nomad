import pytest

from nomad.metainfo.util import MDefNotFound, MDefWithoutMetainfo, resolve_m_def
from tests.metainfo.test_metainfo import SectionWithBoth


@pytest.mark.parametrize(
    'identifier, expected',
    [
        pytest.param(
            'tests.metainfo.test_metainfo.SectionWithBoth',
            SectionWithBoth.m_def,
            id='section',
        ),
        pytest.param(
            'tests.metainfo.test_metainfo.SectionWithBoth.quantity',
            SectionWithBoth.quantity,
            id='quantity',
        ),
        pytest.param(
            'tests.metainfo.test_metainfo.SectionWithBoth.subsection_norepeat',
            SectionWithBoth.subsection_norepeat,
            id='subsection',
        ),
    ],
)
def test_resolve_m_def(identifier, expected):
    assert resolve_m_def(identifier) == expected


@pytest.mark.parametrize(
    'identifier, expected_error',
    [
        pytest.param('nonexistent', MDefNotFound, id='non-existent-module'),
        pytest.param('non.existent', MDefNotFound, id='non-existent-class'),
        pytest.param(
            'tests.metainfo.test_metainfo.SectionWithBoth.nonexistent',
            MDefNotFound,
            id='non-existent-property',
        ),
        pytest.param(
            'nomad.metainfo.data_type.Datatype',
            MDefWithoutMetainfo,
            id='valid-class-no-m_def',
        ),
    ],
)
def test_resolve_m_def_errors(identifier, expected_error):
    with pytest.raises(expected_error):
        resolve_m_def(identifier)
