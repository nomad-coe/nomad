import pytest
import numpy as np
import pint
from nomad.units import ureg
from nomad.parsing.file_parser import (
    TextParser,
    Quantity,
    ParsePattern,
    XMLParser,
    BasicParser,
    FileParser,
)
from nomad.datamodel.metainfo.simulation.system import Atoms
from nomad.datamodel import EntryArchive
from nomad.datamodel.metainfo.simulation.run import Run
from nomad.datamodel.metainfo.simulation.calculation import Calculation


class TestFileParser:
    @pytest.fixture(scope='class')
    def calculation_parser(self):
        class Parser(FileParser):
            def parse(self, key):
                self._results = {'time_calculation': 2.0}

        return Parser()

    @pytest.fixture(scope='class')
    def run_parser(self, calculation_parser):
        class Parser(FileParser):
            def parse(self, key):
                self._results = {'clean_end': True, 'calculation': [calculation_parser]}

        return Parser()

    @pytest.fixture(scope='function')
    def text_parser(self):
        return TextParser()

    @pytest.fixture(scope='class')
    def parser(self, run_parser):
        class Parser(FileParser):
            def parse(self, key):
                self._results = {'run': [run_parser]}

        return Parser()

    @pytest.mark.parametrize(
        'mainfile',
        [
            'tests/data/parsers/vasp/vasp.xml',
            'tests/data/parsers/vasp_compressed/vasp.xml.gz',
        ],
    )
    def test_open(self, text_parser, mainfile):
        text_parser.quantities = [
            Quantity('program', r'name="program" type="string">(.+?) *<')
        ]
        text_parser.mainfile = mainfile
        assert text_parser.program == 'vasp'

    def test_get(self, text_parser):
        text_parser.quantities = [
            Quantity(
                'energy',
                r'free +energy +TOTEN += +(\S+) +eV',
                dtype=np.float64,
                repeats=True,
            ),
        ]
        text_parser.mainfile = 'tests/data/parsers/vasp_outcar/OUTCAR'
        assert text_parser.energy[1] == 1.70437998
        assert text_parser.get('energy', unit='eV')[1].magnitude == 1.70437998

    def test_write_to_archive(self, parser):
        for create_section in [True, False]:
            archive = EntryArchive()
            if create_section:
                archive.m_create(Run).m_create(Calculation)
            archive = parser.write_to_archive(archive)
            assert archive.run[0].clean_end
            assert archive.run[0].calculation[0].time_calculation.magnitude == 2.0


class TestTextParser:
    @pytest.fixture(scope='class')
    def mainfile(self):
        return 'tests/data/parsers/exciting/Ag/INFO.OUT'

    @pytest.fixture(scope='class')
    def mainfile2(self):
        return 'tests/data/parsers/exciting/GW/INFO.OUT'

    @pytest.fixture(scope='class')
    def quantity_string(self):
        return dict(
            quantity=Quantity('spin', r'Spin treatment\s*:\s*([\w\-]+)', repeats=False),
            value='spin-unpolarised',
        )

    @pytest.fixture(scope='class')
    def quantity_float(self):
        return dict(
            quantity=Quantity(
                'cell_volume', r'Unit cell volume\s*:\s*([\d\.]+)', repeats=False
            ),
            value=115.0293819379,
        )

    @pytest.fixture(scope='class')
    def quantity_array(self):
        return dict(
            quantity=Quantity(
                'lattice_vector',
                r'Lattice vectors \(cartesian\) :\s*([\d\s\.]+)',
                repeats=False,
                shape=(3, 3),
            ),
            value=np.array(
                [[3.86005, 3.86005, 0], [3.86005, 0, 3.86005], [0, 3.86005, 3.86005]]
            ),
        )

    @pytest.fixture(scope='class')
    def quantity_repeats(self):
        return dict(
            quantity=Quantity(
                'total_energy', r'Total energy\s*:\s*([\d\.\-]+)', repeats=True
            ),
            value=np.array(
                [
                    -5307.34855605,
                    -5313.90710687,
                    -5315.97055490,
                    -5316.38701749,
                    -5317.59994092,
                    -5317.26163104,
                    -5317.26791647,
                    -5317.26750374,
                    -5317.26724651,
                    -5317.26725951,
                    -5317.26726114,
                    -5317.26726119,
                    -5317.26726118,
                ]
            ),
        )

    @pytest.fixture(scope='class')
    def quantity_with_unit(self):
        return dict(
            quantity=Quantity(
                'wall_time',
                r'Wall time \((?P<__unit>\w+)\)\s*:\s*([\d\.]+)',
                repeats=True,
            ),
            value=[
                pint.Quantity(v, 'seconds')
                for v in [
                    3.55,
                    5.32,
                    7.09,
                    8.84,
                    10.58,
                    12.33,
                    14.09,
                    15.84,
                    17.58,
                    19.33,
                    21.09,
                    22.91,
                ]
            ],
        )

    @pytest.fixture(scope='function')
    def parser(self, mainfile):
        return TextParser(mainfile=mainfile)

    def test_mainfile_setter(self, parser, mainfile2):
        parser.quantities = [
            Quantity('time', r'Time \(hh:mm:ss\)\s*:\s*([\d:]+)', repeats=False)
        ]
        assert parser.get('time') == '19:10:23'
        parser.mainfile = mainfile2
        assert parser.get('time') == '08:24:03'

    def test_constructor(self, mainfile, quantity_string):
        class TestParser(TextParser):
            def __init__(self, **kwargs):
                super().__init__(**kwargs)

            def init_quantities(self):
                self.quantities = [quantity_string['quantity']]

        parser = TestParser(mainfile=mainfile)
        assert parser.get(quantity_string['quantity'].name) == quantity_string['value']

    def test_copy(self, parser, quantity_string):
        parser.quantities = [quantity_string['quantity']]
        parser.parse()
        parser2 = parser.copy()
        assert parser2.mainfile == parser.mainfile
        assert parser2.quantities == parser.quantities
        assert parser2._results != parser._results

    def test_findall(self, parser, quantity_string, quantity_float, quantity_repeats):
        parser.quantities = [
            q['quantity'] for q in [quantity_string, quantity_float, quantity_repeats]
        ]
        assert parser.findall
        spin = parser.get(quantity_string['quantity'].name)
        volume = parser.get(quantity_float['quantity'].name)
        energies = parser.get(quantity_repeats['quantity'].name)

        parser_finditer = parser.copy()
        parser_finditer.findall = False
        assert parser_finditer._results is None
        assert parser_finditer.get(quantity_string['quantity'].name) == spin
        assert parser_finditer.get(quantity_float['quantity'].name) == volume
        assert parser_finditer.get(quantity_repeats['quantity'].name) == energies

    def test_finditer(
        self, parser, quantity_string, quantity_float, quantity_with_unit
    ):
        # dynamic parsing
        parser.quantities = [
            q['quantity'] for q in [quantity_string, quantity_float, quantity_with_unit]
        ]
        parser.findall = False
        count = 0
        for q in [quantity_string, quantity_float, quantity_with_unit]:
            count += 1
            assert parser.get(q['quantity'].name) == q['value']
            assert len(parser._results) == count

    def test_quantity_sub_parser(self, parser):
        quantity_species = Quantity(
            'species',
            r'cies :\s*([\s\S]+?)(?:Spe|Total)',
            repeats=True,
            sub_parser=TextParser(
                quantities=[
                    Quantity('name', r'name\s*:\s*(\w+)', repeats=False),
                    Quantity('mass', r'atomic mass\s*:\s*([\d\.]+)', repeats=False),
                ]
            ),
        )

        quantity_initialization = Quantity(
            'initialization',
            r'Starting initialization([\s\S]+?)Ending initialization',
            repeats=False,
            sub_parser=TextParser(
                quantities=[
                    Quantity(
                        'k_point_grid', r'k\-point grid\s*:\s*([\d ]+)', repeats=False
                    ),
                    quantity_species,
                ]
            ),
        )

        quantity_scf = Quantity(
            'scf',
            r'Self\-consistent loop started([\s\S]+?)Self\-consistent loop stopped',
            repeats=True,
            sub_parser=TextParser(
                quantities=[
                    Quantity(
                        'iteration', r'SCF iteration number\s*:\s*(\d+)', repeats=True
                    )
                ]
            ),
        )

        parser.quantities = [
            quantity_initialization,
            quantity_scf,
            Quantity(
                'total_time',
                r'Total time spent \(seconds\)\s*:\s*([\d.]+)',
                repeats=False,
            ),
        ]

        initialization = parser.get('initialization')

        assert (initialization.get('k_point_grid') == np.array([4, 4, 4])).all()
        species = initialization.get('species')
        assert len(species) == 1
        assert species[0].get('name') == 'silver'
        assert species[0].get('mass') == 196631.6997

        scf = parser.get('scf')
        assert len(scf) == 1
        assert len(scf[0].get('iteration')) == 12

        assert parser.get('total_time') == 22.4

    def test_block_short(self, parser, quantity_repeats):
        parser.quantities = [
            Quantity(
                'scf',
                r'SCF iteration number\s*:\s*\d+([\s\S]+?)Wall time',
                repeats=True,
                sub_parser=TextParser(quantities=[quantity_repeats.get('quantity')]),
            )
        ]

        scf = parser.get('scf')
        assert len(scf) == 12
        energies = quantity_repeats.get('value')
        # total_energy repeats is deliberately set to True to confirm that
        # only one energy per scf block is read
        for i in range(len(scf)):
            assert scf[i].get(quantity_repeats.get('quantity').name) == [energies[i]]

    def test_get_default(self, parser, quantity_float):
        parser.quantities = [quantity_float['quantity']]
        volume = parser.get('volume', 10.00)
        assert volume == 10.00

    def test_get_unit(self, parser, quantity_float):
        parser.quantities = [quantity_float['quantity']]
        volume = parser.get('cell_volume', unit='angstrom')
        assert isinstance(volume, pint.Quantity)
        assert volume.units == 'angstrom'
        assert volume.magnitude == quantity_float['value']

    def test_quantity_unit(
        self, parser, quantity_with_unit, quantity_repeats, quantity_string
    ):
        parser.quantities = [
            q['quantity']
            for q in [quantity_with_unit, quantity_repeats, quantity_string]
        ]
        for q in [quantity_with_unit, quantity_string, quantity_repeats]:
            equal = parser.get(q['quantity'].name) == q['value']
            if isinstance(equal, np.ndarray):
                equal = equal.all()
            assert equal

    def test_quantity_conversion(self, parser, quantity_float):
        quantity = quantity_float.get('quantity')
        quantity.convert = False
        parser.quantities = [quantity]
        assert parser.get(quantity.name) == str(quantity_float.get('value'))

        quantity.convert = True
        parser = parser.copy()
        parser.quantities = [quantity]
        assert parser.get(quantity.name) == quantity_float.get('value')

    def test_quantity_parse_pattern(self, parser):
        parser.quantities = [
            Quantity(
                'BZ_volume',
                ParsePattern(key='Brillouin zone volume', value='re_float'),
                repeats=False,
            ),
            Quantity(
                'n_crystal_symmetry',
                ParsePattern(key='Number of crystal symmetries', value='re_int'),
                repeats=False,
            ),
            Quantity(
                'g_vector_size',
                ParsePattern(key='G-vector grid sizes', value='re_int_array'),
                repeats=False,
            ),
            Quantity(
                'smearing',
                ParsePattern(key='Smearing scheme', value='re_string'),
                repeats=False,
            ),
        ]

        assert parser.get('BZ_volume') == 2.1564074262
        assert parser.get('n_crystal_symmetry') == 48
        assert (parser.get('g_vector_size') == np.array([24, 24, 24])).all()
        assert parser.get('smearing') == 'Gaussian'

    def test_quantity_str_operation(self, parser):
        def str_to_max_lm(string):
            val = [v.split(':') for v in string.split('\n') if v]
            res = {v[0].strip(): v[1] for v in val}
            return res

        parser.quantities = [
            Quantity(
                'max_lm',
                r'Maximum angular momentum used for([\s\S]+?inner part of muffin\-tin\s*:\s*\d+)',
                repeats=False,
                str_operation=str_to_max_lm,
            )
        ]

        max_lm = parser.get('max_lm')
        assert max_lm.get('APW functions') == 8
        assert max_lm.get('computing H and O matrix elements') == 8
        assert max_lm.get('potential and density') == 8
        assert max_lm.get('inner part of muffin-tin') == 2

    def test_quantity_repeats(self, parser, quantity_repeats, mainfile):
        quantity = quantity_repeats.get('quantity')
        parser.quantities = [quantity]
        assert isinstance(parser.get(quantity.name), list)
        quantity.repeats = False
        parser.mainfile = mainfile
        parser.quantities = [quantity]
        assert isinstance(parser.get(quantity.name), float)

    def test_quantity_unit_array(self, parser, quantity_array):
        quantity = quantity_array.get('quantity')
        quantity.unit = 'angstrom'
        parser.quantities = [quantity]
        lattice_vector = parser.get(quantity.name)
        assert isinstance(lattice_vector, pint.Quantity)
        assert lattice_vector.units == 'angstrom'
        assert (lattice_vector.magnitude == quantity_array.get('value')).all()

    def test_quantity_metainfo(self, parser):
        quantity = Quantity(
            Atoms.lattice_vectors,
            r'Lattice vectors \(cartesian\) :\s*([\d\s\.]+)',
            repeats=False,
        )

        parser.quantities = [quantity]
        lattice_vectors = parser.get(quantity.name)
        assert list(lattice_vectors.shape) == Atoms.lattice_vectors.shape
        assert lattice_vectors.dtype == Atoms.lattice_vectors.type


class TestXMLParser:
    @pytest.fixture(scope='class')
    def mainfile(self):
        return 'tests/data/parsers/vasp/vasp.xml'

    @pytest.fixture(scope='function')
    def parser(self, mainfile):
        return XMLParser(mainfile)

    def test_constructor(self, mainfile):
        class TestParser(XMLParser):
            def __init__(self, **kwargs):
                super().__init__(**kwargs)

        test_parser = TestParser(mainfile=mainfile)

        incar = dict(zip(test_parser.get('incar/i/name'), test_parser.get('incar/i')))
        assert incar['SYSTEM'] == 'SrTiO3'
        assert incar['ISMEAR'] == 0
        assert incar['SIGMA'] == 0.1

    def test_parse(self, parser):
        k_points = parser.get('kpoints/varray[1]/v')
        assert isinstance(k_points, np.ndarray)
        assert k_points.shape == (35, 3)

        weights = parser.get('kpoints/varray[2]/v')
        assert isinstance(weights, np.ndarray)
        assert weights.shape == (35,)

        sc_energies = parser.get('calculation/scstep/energy/i')
        assert len(sc_energies) == 13

    def test_conversion(self, parser, mainfile):
        parser.convert = False

        assert parser.get('atominfo/atoms') == ['       5 ']
        assert parser.get('structure[1]/crystal[1]/varray[1]/v[1]') == [
            '       4.00419668       0.00000000       0.00000000 '
        ]

        parser.mainfile = mainfile
        assert parser.get('atominfo/atoms', convert=True) == 5
        parser.convert = True
        assert (
            parser.get('structure[1]/crystal[1]/varray[1]/v[1]')
            == np.array([4.00419668, 0.0, 0.0])
        ).all()


class TestBasicParser:
    @pytest.fixture(scope='class')
    def onetep_parser(self):
        re_f = r'\-*\d+\.\d+E*\-*\+*\d+'
        return BasicParser(
            'ONETEP',
            units_mapping=dict(energy=ureg.hartree, length=ureg.bohr),
            auxilliary_files=r'([\w\-]+\.dat)',
            program_version=r'Version\s*([\d\.]+)',
            lattice_vectors=r'\%block lattice_cart\s*([\s\S]+?)\%endblock lattice_cart',
            atom_labels_atom_positions=rf'\%block positions\_abs\s*(\w+\s+{re_f}\s+{re_f}\s+{re_f}[\s\S]+?)\%endblock positions\_abs',
            XC_functional=r'xc\_functional\s*\:\s*(\w+)',
            energy_total=rf'Total energy\s*=\s*({re_f})\s*Eh',
        )

    def test_onetep_parser(self, onetep_parser):
        archive = EntryArchive()
        onetep_parser.parse(
            'tests/data/parsers/onetep/fluor/12-difluoroethane.out', archive, None
        )

        assert archive.run[0].program.version == '4.5.3.32'
        assert len(archive.run[0].calculation) == 4
        sec_system = archive.run[0].system[0]
        assert sec_system.atoms.labels[7] == 'H'
        assert np.shape(sec_system.atoms.positions) == (8, 3)
