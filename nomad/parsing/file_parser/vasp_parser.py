import numpy as np

import runschema
from nomad.parsing.file_parser.mapping_parser import MappingAnnotationModel, XMLParser


class Program(runschema.run.Program):
    version = runschema.run.Program.version
    version.m_annotations = dict(
        xml=MappingAnnotationModel(
            # path=[".i[?_name=='version']"],
            operator=(
                'get_version',
                [
                    ".i[?_name=='version']",
                    ".i[?_name=='subversion']",
                    ".i[?_name=='platform']",
                ],
            ),
        )
    )

    # compilation_datetime = runschema.run.Program.compilation_datetime
    # compilation_datetime.m_annotations = dict(
    #     xml=XMLAnnotation(operator=[('get_compilation_datetime', ["i[?_name=='']"])])
    # )


class BandEnergies(runschema.calculation.BandEnergies):
    n_spin_channels = runschema.calculation.BandEnergies.n_spin_channels
    n_spin_channels.m_annotations = dict(
        xml=MappingAnnotationModel(path='length(.array.set.set)')
    )

    n_kpoints = runschema.calculation.BandEnergies.n_kpoints
    n_kpoints.m_annotations = dict(
        xml=MappingAnnotationModel(path='length(.array.set.set[0].set)')
    )

    energies = runschema.calculation.BandEnergies.energies
    energies.m_annotations = dict(
        xml=MappingAnnotationModel(
            operator=(
                'get_eigenvalues_energies',
                [
                    '.array.set.set[].set[].r',
                    'length(.array.set.set)',
                    'length(.array.set.set[0].set)',
                ],
            ),
        )
    )


class Calculation(runschema.calculation.Calculation):
    eigenvalues = runschema.calculation.Calculation.eigenvalues
    eigenvalues.m_annotations = dict(xml=MappingAnnotationModel(path='.eigenvalues'))


class Run(runschema.run.Run):
    program = runschema.run.Run.program
    program.m_annotations = dict(xml=MappingAnnotationModel(path='.generator'))

    calculation = runschema.run.Run.calculation
    calculation.m_annotations = dict(xml=MappingAnnotationModel(path='.calculation'))


runschema.run.Run.m_def.m_annotations = dict(
    xml=MappingAnnotationModel(path='modeling')
)


class VASPXMLParser(XMLParser):
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
