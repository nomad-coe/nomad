from copy import deepcopy
from typing import Any

import numpy as np
import pytest

from nomad.datamodel import ArchiveSection
from nomad.datamodel.metainfo.annotations import Mapper as MapperAnnotation
from nomad.metainfo import Quantity, SubSection
from nomad.parsing.file_parser.mapping_parser import (
    MAPPING_ANNOTATION_KEY,
    Data,
    HDF5Parser,
    Mapper,
    MappingParser,
    MetainfoParser,
    Path,
    PathParser,
    TextParser,
    Transformer,
    XMLParser,
)
from nomad.parsing.file_parser.text_parser import Quantity as TextQuantity
from nomad.parsing.file_parser.text_parser import TextParser as TextFileParser


class BSection(ArchiveSection):
    v = Quantity(type=np.float64, shape=[2, 2])
    v.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        xml=MapperAnnotation(mapper='.v'),
        hdf5=MapperAnnotation(mapper=('get_v', ['.v[0].d'])),
    )

    v2 = Quantity(type=str)
    v2.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        xml=MapperAnnotation(mapper='.c[0].d[1]'),
        hdf5=MapperAnnotation(mapper='g.v[-2]'),
        text=MapperAnnotation(mapper='version'),
    )

    v3 = Quantity(type=float)
    v3.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        text=MapperAnnotation(mapper='.energy')
    )


class CSection(ArchiveSection):
    i = Quantity(type=int)
    i.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        xml=MapperAnnotation(mapper='.d'),
        hdf5=MapperAnnotation(mapper='.i | [1]'),
    )

    e = Quantity(type=str)
    e.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        xml=MapperAnnotation(mapper='a.b2.c.e[?"@name"==\'item2\'].k[1] | [0]'),
        hdf5=MapperAnnotation(mapper=('to_string', ['.f[?"@index">=`1`].__value'])),
    )

    g = Quantity(type=np.float64, shape=[2, 5])
    g.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        xml=MapperAnnotation(mapper=('slice', ['a.b2.c.f.g.i'])),
        hdf5=MapperAnnotation(mapper='g.g.c1.d[:2].e[-5:]'),
    )


class B2Section(ArchiveSection):
    c = SubSection(sub_section=CSection)
    c.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        xml=MapperAnnotation(mapper='.c'),
        hdf5=MapperAnnotation(mapper='g.g.c1'),
    )

    b = SubSection(sub_section=BSection)
    b.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        xml=MapperAnnotation(mapper='a.b'),
        hdf5=MapperAnnotation(mapper='.c'),
    )


class ExampleSection(ArchiveSection):
    b = SubSection(sub_section=BSection, repeats=True)
    b.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        xml=MapperAnnotation(mapper='a.b1'),
        hdf5=MapperAnnotation(mapper='.g1'),
        text=MapperAnnotation(mapper='calculation'),
    )
    b2 = SubSection(sub_section=B2Section)
    b2.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        xml=MapperAnnotation(mapper='a.b2'),
        hdf5=MapperAnnotation(mapper='g.g'),
    )


ExampleSection.m_def.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    xml=MapperAnnotation(mapper='a'),
    hdf5=MapperAnnotation(mapper='g'),
    text=MapperAnnotation(),
)


class ExampleXMLParser(XMLParser):
    @staticmethod
    def get_eigenvalues_energies(value, n_spin, n_kpoints):
        array = np.transpose(value)[0].T
        return np.reshape(array, (n_spin, n_kpoints, len(array[0])))

    @staticmethod
    def get_version(version, sub_version, platform):
        return ' '.join([' '.join(s.split()) for s in [version, sub_version, platform]])

    @staticmethod
    def slice(value):
        return np.array(value)[2:]


class ExampleHDF5Parser(HDF5Parser):
    @staticmethod
    def get_v(value):
        return np.array(value)[1:, :2]

    @staticmethod
    def to_string(value):
        return '-'.join([str(n) for n in value])


class ExampleParser(MappingParser):
    def from_dict(self, dct: dict[str, Any]):
        return super().from_dict(dct)  # type: ignore

    def load_file(self) -> Any:
        return super().load_file()

    def to_dict(self, **kwargs) -> dict[str | int, Any]:
        return super().to_dict(**kwargs)

    def slice(self, value):
        return value[1]


@pytest.fixture(scope='module')
def text_parser() -> TextParser:
    outcar_parser = TextFileParser(
        quantities=[
            TextQuantity('version', r'vasp\.([\S]+)'),
            TextQuantity(
                'calculation',
                r'(FREE ENERGIE OF THE ION\-ELECTRON SYSTEM[\s\S]+?entropy.+)',
                repeats=True,
                sub_parser=TextFileParser(
                    quantities=[
                        TextQuantity(
                            'energy',
                            r'free\s*energy\s*TOTEN\s*=\s*([\-\d\.]+)',
                            dtype=float,
                        )
                    ]
                ),
            ),
        ]
    )
    return TextParser(
        text_parser=outcar_parser, filepath='tests/data/parsers/vasp_outcar/OUTCAR'
    )


@pytest.fixture(scope='module')
def xml_parser() -> ExampleXMLParser:
    return ExampleXMLParser(filepath='tests/data/parsing/file_parser/test.xml')


@pytest.fixture(scope='module')
def hdf5_parser() -> ExampleHDF5Parser:
    return ExampleHDF5Parser(filepath='tests/data/parsing/file_parser/test.h5')


@pytest.fixture(scope='module')
def archive_parser() -> MetainfoParser:
    return MetainfoParser()


@pytest.fixture(scope='module')
def data():
    return {
        'a': {
            'b': [
                {'c': {'d': 1, 'e': 'x'}, 'f': [1.0, 2.0]},
                {'c': 2, 'd': [{'e': 'y', 'f': np.eye(2)}]},
            ],
            'c': [
                {'n': 'x', 'v': 1},
                {'n': 'y', 'v': 2},
            ],
        },
        'b': [
            {
                'c': [
                    {'d': 3, 'e': [[1, 2], [3, 4]]},
                    {'d': 4, 'e': [[1, 0], [2, 0]]},
                ]
            },
            {'c': {'d': 1, 'e': 'z'}},
        ],
        'c': [
            {'n': 1, 'v': 'a'},
            {'n': 2, 'v': 'b'},
        ],
    }


def assert_equal(v1, v2):
    if isinstance(v1, dict):
        for key, val in v1.items():
            assert key in v2
            assert_equal(val, v2[key])
    elif isinstance(v1, list):
        assert isinstance(v2, list)
        for n, v in enumerate(v1):
            assert_equal(v, v2[n])
    else:
        equal = v1 == v2
        assert equal.all() if isinstance(equal, np.ndarray) else equal


class TestPath:
    @pytest.mark.parametrize(
        'path, result',
        [
            pytest.param('a.b[1].c', 2),
            pytest.param('b[0].c[0:2].d', [3, 4]),
            pytest.param('b[0].c[0:2].e', [[[1, 2], [3, 4]], [[1, 0], [2, 0]]]),
            pytest.param('a.b[1].d[0].f', np.eye(2)),
            pytest.param('b[1].c.e', 'z'),
            pytest.param('a[1].b.c.d', None),
            pytest.param("a.c[?n=='x'].v | [0]", 1),
            pytest.param('c[?n>=`1`].v | [1]', 'b'),
        ],
    )
    def test_get_data(self, path, result, data):
        path = Path(path=path)
        value = path.get_data(data)
        assert_equal(value, result)

    @pytest.mark.parametrize(
        'path, data, target, result',
        [
            pytest.param('a.b', 'x', {}, 'x'),
            pytest.param('a[1:3].b.c[0].d[0:3]', 'y', {}, [['y'] * 3] * 2),
            pytest.param('a.b[0:3].c[1:2]', [[1], [2], [3]], {}, [[1], [2], [3]]),
            pytest.param('a[1:4].b', ['a', 'b', 'c'], {}, ['a', 'b', 'c']),
            pytest.param('a[0:2].b.c', [1, 2, 3], {}, [[1, 2, 3]] * 2),
            pytest.param(
                'a[0].b[0:2].c[0]',
                [['x'], ['y']],
                {'a': [{'b': [{'c': ['a']}, {'c': ['b']}]}]},
                ['x', 'y'],
            ),
            pytest.param(
                'a[0].b[0:2]',
                [['x'], ['y']],
                {'a': [{'b': ['a', 'b']}]},
                [['x'], ['y']],
            ),
        ],
    )
    def test_set_data(self, path, data, target, result):
        path = Path(path=path)
        path.set_data(data, target)
        value = path.get_data(target)
        assert_equal(value, result)

    def test_set_data_nested_update_mode_overrides_parent(self):
        path = Path(path='items')
        target = {
            'items': [
                {
                    'entries': [
                        {
                            'kind': 'existing',
                            'value': 1,
                        }
                    ]
                }
            ]
        }
        data = [
            {
                '.entries': [
                    {
                        '.kind': 'incoming',
                        '.label': 'x',
                    }
                ]
            }
        ]

        path.set_data(
            data,
            target,
            update_mode={
                '__update_mode': 'merge',
                '.entries': {'__update_mode': 'append'},
            },
        )

        item_data = target['items'][0]
        entries = item_data.get('entries', item_data.get('.entries', []))
        assert len(entries) == 2
        assert any('label' in section or '.label' in section for section in entries)
        assert any('value' in section or '.value' in section for section in entries)
        assert {
            section.get('kind', section.get('.kind'))
            for section in entries
            if 'kind' in section or '.kind' in section
        } == {
            'existing',
            'incoming',
        }

    def test_set_data_merge_last_scalar_list(self):
        path = Path(path='a.b')
        target = {'a': {'b': [10, 20]}}

        path.set_data([1, 2, 3, 4, 5], target, update_mode='merge@last')

        assert target['a']['b'] == [1, 2, 3, 10, 20]

    def test_set_data_merge_last_scalar_list_current_longer(self):
        path = Path(path='a.b')
        target = {'a': {'b': [10, 20, 30, 40, 50]}}

        path.set_data([1, 2], target, update_mode='merge@last')

        assert target['a']['b'] == [10, 20, 30, 40, 50]

    def test_set_data_merge_invalid_index_raises(self):
        path = Path(path='a.b')
        target = {'a': {'b': [10, 20]}}

        with pytest.raises(ValueError, match='merge index must be an integer'):
            path.set_data([1, 2, 3], target, update_mode='merge@foo')


class TestMapper:
    @pytest.mark.parametrize(
        'dct, expected',
        [
            pytest.param(
                dict(source='a', target='b', mapper='v'),
                Transformer(
                    source=Data(path=Path(path='a')),
                    target=Data(path=Path(path='b')),
                    function_args=[Path(path='v')],
                ),
            ),
            pytest.param(
                dict(source='a', path='.v', path_parser='jsonpath_ng'),
                Transformer(
                    source=Data(
                        path=Path(
                            path='a', parser=PathParser(parser_name='jsonpath_ng')
                        ),
                        path_parser=PathParser(parser_name='jsonpath_ng'),
                    ),
                    function_args=[
                        Path(
                            parser=PathParser(parser_name='jsonpath_ng'),
                            path='.v',
                            parent=Path(
                                path='a', parser=PathParser(parser_name='jsonpath_ng')
                            ),
                        )
                    ],
                ),
            ),
            pytest.param(
                dict(target='b', source='a', mapper=('eval', ['.a', 'b'])),
                Transformer(
                    source=Data(path=Path(path='a')),
                    target=Data(path=Path(path='b')),
                    function_name='eval',
                    function_args=[
                        Path(path='.a', parent=Path(path='a')),
                        Path(path='b'),
                    ],
                ),
            ),
            pytest.param(
                dict(
                    source=Data(
                        transformer=Transformer(function_args=[Path(path='a')])
                    ),
                    function_name='eval',
                    function_args=[Path(path='.b')],
                ),
                Transformer(
                    source=Data(
                        transformer=Transformer(function_args=[Path(path='a')])
                    ),
                    function_name='eval',
                    function_args=[Path(path='.b', parent=Path(path='a'))],
                ),
            ),
            pytest.param(
                dict(
                    source='a',
                    remove=True,
                    mapper=[
                        dict(path='.c', source='.b'),
                        dict(
                            path_parser='jsonpath_ng',
                            mapper=['eval', ['.x']],
                            remove=False,
                        ),
                        dict(
                            path_parser='jsonpath_ng',
                            mapper=[dict(mapper='.d', remove=True)],
                            source='.b',
                            remove=False,
                        ),
                    ],
                ),
                Mapper(
                    source=Data(path=Path(path='a')),
                    remove=True,
                    mappers=[
                        Transformer(
                            function_args=[
                                Path(
                                    path='.c',
                                    parent=Path(path='.b', parent=Path(path='a')),
                                )
                            ],
                            source=Data(
                                path=Path(path='.b', parent=Path(path='a')),
                                parent=Path(path='a'),
                            ),
                            remove=True,
                        ),
                        Transformer(
                            function_name='eval',
                            function_args=[
                                Path(
                                    path='.x',
                                    parser=PathParser(parser_name='jsonpath_ng'),
                                    parent=Path(path='a'),
                                )
                            ],
                            remove=False,
                        ),
                        Mapper(
                            mappers=[
                                Transformer(
                                    function_args=[
                                        Path(
                                            path='.d',
                                            parent=Path(
                                                path='.b',
                                                parser=PathParser(
                                                    parser_name='jsonpath_ng'
                                                ),
                                                parent=Path(path='a'),
                                            ),
                                        )
                                    ],
                                    remove=True,
                                )
                            ],
                            source=Data(
                                path=Path(
                                    path='.b',
                                    parser=PathParser(parser_name='jsonpath_ng'),
                                    parent=Path(path='a'),
                                ),
                                path_parser=PathParser(parser_name='jsonpath_ng'),
                                parent=Path(path='a'),
                            ),
                            remove=False,
                        ),
                    ],
                ),
            ),
        ],
    )
    def test_from_dict(self, dct, expected):
        def assert_mappers_equal(m1, m2):
            assert isinstance(m1, type(m2))
            assert m1.source == m2.source
            assert m1.target == m2.target
            assert m1.remove == m2.remove
            assert m1.indices == m2.indices
            if isinstance(m1, Mapper):
                for n, sm1 in enumerate(m1.mappers):
                    assert_mappers_equal(sm1, m2.mappers[n])
            elif isinstance(m1, Transformer):
                assert m1.function_name == m2.function_name
                for n, arg in enumerate(m1.function_args):
                    assert arg == m2.function_args[n]

        mapper = Mapper.from_dict(dct)
        assert_mappers_equal(mapper, expected)

    @pytest.mark.parametrize('remove', [True, False])
    @pytest.mark.parametrize(
        'mapper, expected',
        [
            pytest.param(
                Mapper(
                    mappers=[
                        Transformer(
                            source=Data(
                                path=Path(
                                    path='a.b',
                                    parser=PathParser(parser_name='jsonpath_ng'),
                                )
                            ),
                            function_args=[
                                Path(
                                    path='.f',
                                    parser=PathParser(parser_name='jsonpath_ng'),
                                )
                            ],
                            target=Data(path=Path(path='x')),
                        )
                    ]
                ),
                dict(x=[1.0, 2.0]),
            ),
            pytest.param(
                Mapper(
                    mappers=[
                        Transformer(
                            source=Data(
                                transformer=Transformer(
                                    function_args=[Path(path='a.b')]
                                )
                            ),
                            function_args=[Path(path='.f')],
                            target=Data(path=Path(path='x')),
                        )
                    ],
                ),
                dict(x=[1.0, 2.0]),
            ),
            pytest.param(
                Mapper(
                    mappers=[
                        Transformer(
                            source=Data(
                                transformer=Transformer(
                                    function_args=[Path(path='b[0].c')]
                                )
                            ),
                            function_args=[Path(path='.e')],
                            function_name='slice',
                            target=Data(path=Path(path='x')),
                        )
                    ]
                ),
                dict(x=[3.0, 4.0]),
            ),
            pytest.param(
                Mapper(
                    mappers=[
                        Mapper(
                            source=Data(path=Path(path='a')),
                            mappers=[
                                Mapper(
                                    mappers=[
                                        Transformer(
                                            function_args=[Path(path='.d[0].e')],
                                            target=Data(path=Path(path='z')),
                                        )
                                    ],
                                    source=Data(path=Path(path='.b')),
                                    indices=None,
                                    target=Data(path=Path(path='y')),
                                )
                            ],
                            target=Data(path=Path(path='x')),
                        ),
                        Transformer(
                            function_args=[Path(path='c[?n==`2`].v | [0]')],
                            target=Data(path=Path(path='x2')),
                        ),
                    ]
                ),
                dict(x=dict(y=dict(z='y')), x2='b'),
            ),
        ],
    )
    def test_get_data(self, data, remove, mapper, expected):
        source = deepcopy(data)
        parser = ExampleParser(data=source)
        mapper.remove = remove
        result = mapper.get_data(source, parser)
        assert_equal(expected, result)
        if remove:
            assert not mapper.get_data(source, parser)


class TestMappingParser:
    def test_set_data_nested_update_mode_per_key(self, monkeypatch):
        parser = ExampleParser(data={})
        update_modes: list[tuple[str, str | None]] = []

        def fake_set_data(self, data, target, **kwargs):
            update_modes.append((self.path, kwargs.get('update_mode')))
            return data

        monkeypatch.setattr(Path, 'set_data', fake_set_data)
        parser.set_data(
            {'a': {'b': 1, 'c': 2}},
            {},
            update_mode={
                '__update_mode': 'merge',
                'a': {
                    '__update_mode': 'merge',
                    '.b': {'__update_mode': 'replace'},
                    '.c': {'__update_mode': 'append'},
                },
            },
        )

        mode_map = {path: mode for path, mode in update_modes}
        assert mode_map['a']['__update_mode'] == 'merge'
        assert mode_map['b']['__update_mode'] == 'replace'
        assert mode_map['c']['__update_mode'] == 'append'

    def test_convert_xml_to_archive(self, xml_parser, archive_parser):
        archive_parser.annotation_key = 'xml'
        archive_parser.data_object = ExampleSection(b=[BSection(v=np.eye(2))])

        xml_parser.convert(archive_parser, update_mode='append')
        archive = archive_parser.data_object
        assert len(archive.b) == 3
        assert archive.b[0].v[0][0] == 1.0
        assert archive.b[1].v[1][0] == 3.0
        assert archive.b[2].v[1][1] == 8.0
        assert archive.b[2].v2 == 'b'
        assert archive.b2.c.i == 1
        assert archive.b2.c.e == 'f4'
        assert archive.b2.c.g[1][2] == 8
        xml_parser.close()

    def test_convert_archive_to_xml(self, xml_parser, archive_parser):
        archive_parser.data_object = ExampleSection(b=[BSection(v=np.eye(2))])
        xml_parser.mapper = Mapper(
            mappers=[
                Mapper(
                    target=Data(path=Path(path='a')),
                    mappers=[
                        Mapper(
                            target=Data(path=Path(path='.b1')),
                            mappers=[
                                Transformer(
                                    function_args=[Path(path='.v')],
                                    target=Data(path=Path(path='.v')),
                                )
                            ],
                            source=Data(path=Path(path='b')),
                        )
                    ],
                )
            ],
        )
        xml_parser.filepath = None
        archive_parser.convert(xml_parser)
        assert xml_parser.data_object.findall('b1')[0].findall('v')[1].text == '0.0 1.0'
        xml_parser.close()

    def test_convert_hdf5_to_archive(self, hdf5_parser, archive_parser):
        archive_parser.annotation_key = 'hdf5'
        archive_parser.data_object = ExampleSection(b=[BSection(v=np.eye(2))])
        hdf5_parser.convert(archive_parser, update_mode='merge')
        archive = archive_parser.data_object
        assert archive.b[0].v[1][1] == 1.0
        assert archive.b[0].v2 == 'y'
        assert archive.b2.c.i == 6
        assert archive.b2.c.e == '2-1'
        assert archive.b2.c.g[1][3] == 9
        assert archive.b2.b.v[0][1] == 1
        hdf5_parser.close()

    def test_convert_text_to_archive(self, text_parser, archive_parser):
        archive_parser.annotation_key = 'text'
        archive_parser.data_object = ExampleSection(b=[BSection(v=np.eye(2))])
        text_parser.convert(archive_parser, update_mode='replace')
        archive = archive_parser.data_object
        assert len(archive.b) == 3
        assert archive.b[0].v2 == '5.3.2'
        assert archive.b[2].v3 == -7.14173545
        text_parser.close()

    def test_from_dict_polymorphic_subsections(self):
        """Test that from_dict resolves types per-item for polymorphic subsections."""
        from nomad.metainfo import MSection, Quantity, SubSection

        # Define base section and two concrete types
        class BaseItem(MSection):
            name = Quantity(type=str)

        class ItemTypeA(BaseItem):
            property_a = Quantity(type=str)

        class ItemTypeB(BaseItem):
            property_b = Quantity(type=str)

        class Container(MSection):
            items = SubSection(sub_section=BaseItem, repeats=True)

        # Create parser with Container
        parser = MetainfoParser()
        parser.data_object = Container()

        # Get qualified names from the section definitions
        type_a_qname = ItemTypeA.m_def.qualified_name()
        type_b_qname = ItemTypeB.m_def.qualified_name()

        # Test dict with heterogeneous list: both type A and type B items
        data = {
            'items': [
                {
                    'm_def': type_a_qname,
                    'name': 'first',
                    'property_a': 'value_a',
                },
                {'m_def': type_b_qname, 'name': 'second', 'property_b': 'value_b'},
            ]
        }

        # Call from_dict
        parser.from_dict(data)

        # Verify both items were created with correct types
        assert len(parser.data_object.items) == 2
        assert isinstance(parser.data_object.items[0], ItemTypeA)
        assert parser.data_object.items[0].name == 'first'
        assert parser.data_object.items[0].property_a == 'value_a'

        assert isinstance(parser.data_object.items[1], ItemTypeB)
        assert parser.data_object.items[1].name == 'second'
        assert parser.data_object.items[1].property_b == 'value_b'
