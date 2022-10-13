import yaml
import json
import os.path

from nomad.metainfo import Package
from nomad.units import ureg

from tests.parsing.test_parsing import run_singular_parser


def _file(path):
    return os.path.join(
        os.path.dirname(__file__),
        f'../../examples/docs/{path}')


def _load_yaml(path):
    with open(_file(path), 'rt') as f:
        return yaml.safe_load(f)


def _parse_archive(path):
    return run_singular_parser('parsers/archive', _file(path))


def test_python_schema():
    yaml_data = _load_yaml('basic_schema/data.archive.yaml')['data']
    del(yaml_data['m_def'])

    from examples.docs.basic_schema.schema import Sample
    sample = Sample.m_from_dict(yaml_data)
    assert json.dumps(sample.m_to_dict()) == json.dumps(yaml_data)


def test_yaml_schema():
    yaml_package = _load_yaml('basic_schema/schema.archive.yaml')['definitions']
    yaml_data = _load_yaml('basic_schema/data.archive.yaml')['data']
    del(yaml_data['m_def'])

    package = Package.m_from_dict(yaml_package)
    package.init_metainfo()
    composition_def = package.all_definitions['Composition']
    composition = composition_def.section_cls.m_from_dict(yaml_data)
    assert json.dumps(composition.m_to_dict()) == json.dumps(yaml_data)


def test_yaml_data():
    archive = _parse_archive('basic_schema/data.archive.yaml')
    assert archive.data.composition == 'H2O'
    assert archive.data.elements[0].label == 'H'


def test_references():
    archive = _parse_archive('references/periodic_table.archive.yaml')
    assert archive.data.elements[0].label == 'H'

    archive = _parse_archive('references/composition.archive.yaml')
    assert archive.data.elements[0].label == 'H'

    archive = _parse_archive('references/single.archive.yaml')
    assert archive.data.compositions[0].elements[0].label == 'H'


def test_inheritance():
    archive = _parse_archive('inheritance/basic.archive.yaml')
    assert archive.data.pressure == 100 * ureg('Pa')

    archive = _parse_archive('inheritance/specialized.archive.yaml')
    assert archive.data.processes[0].pressure == 100 * ureg('Pa')
    assert archive.data.processes[1].temperature == 342 * ureg('K')


def test_multiple_files():
    archive = _parse_archive('references/multiple_files/schema.archive.yaml')
    assert len(archive.definitions.sections) == 3

    archive = _parse_archive('references/multiple_files/data-and-schema.archive.yaml')
    assert archive.data.elements[0].label == 'H'
    assert archive.data.elements[1].label == 'O'

    archive = _parse_archive('references/multiple_files/data.archive.yaml')
    assert archive.data.solvent.elements[0].label == 'H'
    assert archive.data.solute.elements[0].label == 'Na'
