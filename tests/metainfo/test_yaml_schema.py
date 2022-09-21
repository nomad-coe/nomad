import numpy as np            # pylint: disable=unused-import
import pytest
import yaml

from nomad.datamodel.data import UserReference, AuthorReference
from nomad.metainfo.metainfo import MTypes
from nomad.utils import strip

from nomad.metainfo import Package, MSection, Quantity, Reference, SubSection, Section, MProxy, MetainfoError

m_package = Package()


class Sample(MSection):
    sample_id = Quantity(
        type=str,
        a_eln=dict(component="StringEditQuantity"),
        description='''
        This is a description with *markup* using [markdown](https://markdown.org).
        It can have multiple lines, because yaml allows to easily do this.
        ''')


class Process(MSection):
    samples = Quantity(type=Reference(Sample.m_def), shape=["*"])
    layers = SubSection(sub_section=Sample.m_def, repeats=True)


class SpecialProcess(Process):
    values = Quantity(type=np.dtype(np.float64), shape=[3, 3])


m_package.__init_metainfo__()


def yaml_to_package(yaml_str):
    yaml_obj = yaml.safe_load(yaml_str)
    package = Package.m_from_dict(yaml_obj)
    return package


yaml_schema_example = strip('''
m_def: 'nomad.metainfo.metainfo.Package'
sections:
  Sample:
    base_section: 'nomad.datamodel.metainfo.measurements.Sample'
    quantities:
      sample_id:
        type: str
        description: |
          This is a description with *markup* using [markdown](https://markdown.org).
          It can have multiple lines, because yaml allows to easily do this.
        m_annotations:
          eln:
            component: StringEditQuantity
  Process:
    quantities:
      samples:
        type: '#/Sample'
        shape: ['*']
    sub_sections:
      layers:
        section_def: '#/Sample'
        repeats: true
  SpecialProcess:
    base_section: '#/Process'
    quantities:
      values:
        type: np.float64
        shape: [3, 3]
''')


def test_yaml_deserialization():
    des_m_package = yaml_to_package(yaml_schema_example)

    sample = m_package['section_definitions'][0]
    process = m_package['section_definitions'][1]
    special_process = m_package['section_definitions'][2]

    des_sample = des_m_package['section_definitions'][0]
    des_process = des_m_package['section_definitions'][1]
    des_special_process = des_m_package['section_definitions'][2]

    def assert_referenced_section(section):
        assert section is not None
        section = section.m_resolved()
        assert isinstance(section, Section)
        assert not isinstance(section, MProxy)

    assert sample.name == des_sample.name == "Sample"

    sample_id = sample['quantities'][0]
    des_sample_id = des_sample['quantities'][0]
    assert sample_id.name == des_sample_id.name == 'sample_id'
    assert sample_id.type == des_sample_id.type == str
    assert sample_id.shape == des_sample_id.shape
    assert sample_id.m_annotations["eln"]["component"] == des_sample_id.m_annotations["eln"]["component"]
    assert sample_id.description == des_sample_id.description.rstrip('\n')

    assert process.name == des_process.name == "Process"

    samples = process['quantities'][0]
    des_samples = des_process['quantities'][0]
    assert samples.name == des_samples.name == 'samples'
    assert samples.shape == des_samples.shape == ["*"]
    assert_referenced_section(des_samples.type.target_section_def)
    assert samples.type.target_section_def.name == des_samples.type.target_section_def.name == 'Sample'

    layers = process['sub_sections'][0]
    des_layers = des_process['sub_sections'][0]
    assert layers.name == des_layers.name == 'layers'
    assert layers.repeats is des_layers.repeats is True
    assert_referenced_section(des_layers.section_def)
    assert layers.section_def.name == des_layers.section_def.name == 'Sample'

    assert special_process.name == des_special_process.name == "SpecialProcess"

    values = special_process['quantities'][0]
    des_values = des_special_process['quantities'][0]
    assert values.name == des_values.name == 'values'
    assert values.shape == des_values.shape == [3, 3]

    base_section = special_process['base_sections'][0]
    des_base_section = des_special_process['base_sections'][0]
    assert_referenced_section(des_base_section)
    assert base_section.name == des_base_section.name == 'Process'

    des_m_package.m_to_dict()


@pytest.mark.parametrize('yaml_schema, expected_error', [
    pytest.param(strip('''
        m_def: 'nomad.metainfo.metainfo.Package'
        sections:
          Process:
            quantities:
              samples:
                type: np.float6
    '''), 'float6 is not a valid numpy type.', id='wrong np type'),
    pytest.param(strip('''
        m_def: 'nomad.metainfo.metainfo.Package'
        sections:
          Process:
            quantities:
              samples:
                type: numpy.int3
    '''), 'int3 is not a valid numpy type.', id='wrong numpy type'),
    pytest.param(strip('''
        m_def: 'nomad.metainfo.metainfo.Package'
        sections:
          Process:
            quantities:
              samples:
                type: float
                m_annotations: eln
    '''), 'The provided m_annotations is of a wrong type. str was provided.', id='wrong m_annotations')
])
def test_errors(yaml_schema, expected_error):
    with pytest.raises(Exception) as exception:
        yaml_to_package(yaml_schema)

    assert isinstance(exception.value, MetainfoError)
    assert exception.value.args[0] == expected_error


def test_sub_section_tree():
    yaml = yaml_to_package('''
      sections:
        Parent:
          sub_sections:
            the_child:
              sub_section:
                quantities:
                  quantity:
                    type: str
    ''')
    reference = yaml_to_package('''
      sections:
        Parent:
          sub_sections:
            the_child:
              section: TheChild
          sections:
            TheChild:
              quantities:
                quantity:
                  type: str
    ''')

    assert yaml.m_to_dict() == reference.m_to_dict()


@pytest.mark.parametrize("eln_type", MTypes.eln.keys())
@pytest.mark.parametrize("eln_component", sum(MTypes.eln_component.values(), []))
def test_datatype_component_annotations(eln_type, eln_component):
    base_schema = '''
              m_def: 'nomad.metainfo.metainfo.Package'
              sections:
                Sample:
                  base_section: 'nomad.datamodel.metainfo.measurements.Sample'
                  quantities:
                    sample_id:
                      type: str
                      m_annotations:
                        eln:
                          component: StringEditQuantity
                Process:
                  quantities:
                    quantity_name:
                      type: quantity_type
                      m_annotations:
                        eln:
                          component: eln_component
            '''

    for quantity_type in MTypes.eln[eln_type]:
        if eln_type == 'reference':
            yaml_schema = base_schema.replace("quantity_type", "'#/Sample'").replace("eln_component", eln_component)
        else:
            yaml_schema = base_schema.replace("quantity_type", quantity_type).replace("eln_component", eln_component)

        if eln_component not in MTypes.eln_component[eln_type]:
            with pytest.raises(Exception) as exception:
                package = yaml_to_package(yaml_schema)
                type_name = quantity_type
                if eln_type == 'number' or eln_type == 'datetime' or eln_type == 'enum' or eln_type == 'reference':
                    process = next(filter(lambda section: section['name'] == 'Process', package['section_definitions']),
                                   None)
                    quantity = process['quantities'][0]
                    if type(quantity.type).__name__ != 'type':
                        type_name = type(quantity.type).__name__
                package.__init_metainfo__()
            assert isinstance(exception.value, MetainfoError)
            assert exception.value.args[0] == 'One constraint was violated: The component `%s` is not compatible with the quantity `%s` of the type `%s`. Accepted components: %s (there are 0 more violations)' \
                % (eln_component, 'quantity_name', type_name, ', '.join(MTypes.eln_component[eln_type]))


yaml_schema_user_author = strip('''
m_def: 'nomad.metainfo.metainfo.Package'
sections:
  Sample:
    base_section: 'nomad.datamodel.metainfo.measurements.Sample'
    quantities:
      my_user:
        type: User
        m_annotations:
          eln:
            component: UserEditQuantity
      my_author:
        type: Author
        m_annotations:
          eln:
            component: AuthorEditQuantity
''')


def test_user_author_yaml_deserialization():
    des_m_package = yaml_to_package(yaml_schema_user_author)
    des_sample = des_m_package['section_definitions'][0]
    des_my_user = des_sample.quantities[0]
    des_my_author = des_sample.quantities[1]

    assert des_my_user.name == 'my_user'
    assert des_my_author.name == 'my_author'
    assert isinstance(des_my_user.type, UserReference)
    assert isinstance(des_my_author.type, AuthorReference)
