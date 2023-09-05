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

from enum import Enum
from typing import Optional, Dict, List
from pydantic import BaseModel, Field

from .models import NomadSettings, Options


_jupyterhub_config_description = '''
    This setting is forwarded to jupyterhub; refer to the jupyterhub
    [documentation](https://jupyterhub.readthedocs.io/en/stable/api/app.html#).
'''


class NORTHToolMaintainer(BaseModel):
    name: str
    email: str


class ReadMode(str, Enum):
    ro = 'ro'
    rw = 'rw'


class NORTHExternalMount(BaseModel):
    host_path: str
    bind: str
    mode: ReadMode = ReadMode.ro


class NORTHTool(BaseModel):
    image: str
    description: str = None
    short_description: str = None
    cmd: str = None
    path_prefix: str = None
    mount_path: str = None
    icon: str = None
    file_extensions: List[str] = []
    maintainer: List[NORTHToolMaintainer] = []
    privileged: bool = False
    external_mounts: List[NORTHExternalMount] = []


class NORTHTools(Options):
    options: Dict[str, NORTHTool] = Field(dict(), description='The available plugin.')


class NORTH(NomadSettings):
    '''
    Settings related to the operation of the NOMAD remote tools hub service *north*.
    '''
    enabled: Optional[str] = Field(description='''
        Enables or disables the NORTH API and UI views. This is independent of
        whether you run a jupyter hub or not.
    ''')
    hub_connect_ip: str = Field(None, description='''
        Overwrites the default hostname that can be used from within a north container
        to reach the host system.

        Typically has to be set for non Linux hosts. Set this to `host.docker.internal`
        on windows/macos.
    ''')
    hub_connect_url: str = Field(None, description=_jupyterhub_config_description)
    hub_ip = Field('0.0.0.0', description=_jupyterhub_config_description)
    docker_network: str = Field(None, description=_jupyterhub_config_description)
    hub_host = Field('localhost', description='''
        The internal host name that NOMAD services use to connect to the jupyterhub API.
    ''')
    hub_port = Field(9000, description='''
        The internal port that NOMAD services use to connect to the jupyterhub API.
    ''')
    jupyterhub_crypt_key: str = Field(None, description=_jupyterhub_config_description)

    nomad_host: str = Field(
        None, description='The NOMAD app host name that spawned containers use.')
    windows = Field(
        True, description='Enable windows OS hacks.')

    tools: NORTHTools = Field(
        NORTHTools(options={
            'jupyter': NORTHTool(
                short_description='Basic jupyter run with an empty notebook or on given notebook file.',
                description='### **Jupyter Notebook**: The Classic Notebook Interface\n\nThe Jupyter Notebook is the original web application for creating and sharing computational documents. It offers a simple, streamlined, document-centric experience.',
                image='gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-remote-tools-hub/jupyterlab:v0.0.1',
                path_prefix='tree',
                mount_path='/home/jovyan',
                icon='jupyter_logo.svg',
                file_extensions=[
                    'ipynb'
                ],
                maintainer=[
                    NORTHToolMaintainer(
                        name='Markus Scheidgen',
                        email='markus.scheidgen@physik.hu-berlin.de'
                    ),
                    NORTHToolMaintainer(
                        name='Some-one Else',
                        email='markus.scheidgen@physik.hu-berlin.de'
                    )
                ]
            ),
            'nionswift': NORTHTool(
                short_description='Run NionSwift to analyze data as well as prepare focus series reconstructions',
                description='Run Nion Swift to analyze data.',
                image='gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-remote-tools-hub/nionswift-webtop',
                privileged=True,
                mount_path='/config',
                maintainer=[
                    NORTHToolMaintainer(
                        name='Sherjeel Shabih',
                        email='sherjeel.shabih@hu-berlin.de'
                    )
                ]
            ),
            'nexustools': NORTHTool(
                description='Includes multiple NeXus tools for visualization and analysis.',
                image='gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-remote-tools-hub/nexus-webtop',
                privileged=True,
                mount_path='/config',
                file_extensions=[
                    'nxs',
                    'nx',
                    'nexus',
                    'hdf5',
                    'hd5',
                    'h5',
                    'hdf',
                    'ipynb'
                ],
                maintainer=[
                    NORTHToolMaintainer(
                        name='Sandor Brockhauser',
                        email='sandor.brockhauser@physik.hu-berlin.de'
                    )
                ]
            ),
            'ellips': NORTHTool(
                short_description='An example for analyzing ellipsometric data.',
                description='This example presents the capabilities of the NOMAD platform to store and standardize ellipsometry data. It shows the generation of a NeXus file according to the [NXellipsometry](https://manual.nexusformat.org/classes/contributed_definitions/NXellipsometry.html#nxellipsometry) application definition and a successive analysis of a SiO2 on Si Psi/Delta measurement.',
                image='gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-remote-tools-hub/ellips-jupyter',
                path_prefix='lab/tree',
                mount_path='/home/jovyan',
                icon='jupyter_logo.svg',
                file_extensions=[
                    'ipynb',
                    'nxs'
                ],
                maintainer=[
                    NORTHToolMaintainer(
                        name='Florian Dobener',
                        email='florian.dobener@physik.hu-berlin.de'
                    ),
                    NORTHToolMaintainer(
                        name='Carola Emminger',
                        email='emminger.carola@physik.hu-berlin.de'
                    )
                ]
            ),
            'mpes': NORTHTool(
                short_description='An example for analyzing mpes data.',
                description='This example presents the capabilities of the NOMAD platform to store and standardize multi photoemission spectroscopy (MPES) experimental data. It contains three major examples:\n\n- Taking a pre-binned file, here stored in a h5 file, and converting it into the standardized MPES NeXus format. There exists a [NeXus application definition for MPES](https://manual.nexusformat.org/classes/contributed_definitions/NXmpes.html#nxmpes) which details the internal structure of such a file.\n- Binning of raw data (see [here](https://www.nature.com/articles/s41597-020-00769-8) for additional resources) into a h5 file and consecutively generating a NeXus file from it.\n- An analysis example using data in the NeXus format and employing the [pyARPES](https://github.com/chstan/arpes) analysis tool to reproduce the main findings of [this paper](https://arxiv.org/pdf/2107.07158.pdf).',
                image='gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-remote-tools-hub/mpes-webtop',
                privileged=True,
                mount_path='/config',
                file_extensions=[
                    'ipynb',
                    'nxs',
                    'h5',
                    'hdf5'
                ],
                maintainer=[
                    NORTHToolMaintainer(
                        name='Florian Dobener',
                        email='florian.dobener@physik.hu-berlin.de'
                    )
                ]
            ),
            'xps': NORTHTool(
                short_description='An example for analyzing XPS data.',
                description='Includes tools for analyzing X-ray Photoelectron Spectroscopy (XPS) spectra and converting SPECS xml files into NeXus.',
                image='gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-remote-tools-hub/xps-jupyter',
                path_prefix='lab/tree',
                icon='jupyter_logo.svg',
                mount_path='/home/jovyan',
                file_extensions=[
                    'ipynb',
                    'nxs',
                    'h5',
                    'hdf5'
                ],
                maintainer=[
                    NORTHToolMaintainer(
                        name='Florian Dobener',
                        email='florian.dobener@physik.hu-berlin.de'
                    ),
                    NORTHToolMaintainer(
                        name='Rubel Mozumder',
                        email='rubel.mozumder@physik.hu-berlin.de'
                    )
                ]
            ),

            'sts': NORTHTool(
                short_description='An example for analyzing SPM (STM /STS) experiment.',
                description="AT this moment, the reader works for two types of experiments, Scanning Tunneling Microscopy (STM) and Scanning Tunneling Spectroscopy (STS) from Scanning Probe Microscopy. It can only transform the data from Nanonis machine generated files into standarised nexus application definition NXsts. The present version of STS reader can handle files from two specific software versions generic 5e and genric 4.5.",
                image='gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-remote-tools-hub/sts-jupyter',
                path_prefix='lab/tree',
                icon='jupyter_logo.svg',
                mount_path='/home/jovyan',
                file_extensions=[
                    'ipynb',
                    'nxs',
                    'h5',
                    'hdf5'
                ],
                maintainer=[
                    NORTHToolMaintainer(
                        name='Rubel Mozumder',
                        email='rubel.mozumder@physik.hu-berlin.de'
                    )
                ]
            ),
            'webtop': NORTHTool(
                description='Baseline webtop image for test',
                image='gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-remote-tools-hub/webtop',
                privileged=True,
                mount_path='/config',
                maintainer=[
                    NORTHToolMaintainer(
                        name='Sherjeel Shabih',
                        email='sherjeel.shabih@hu-berlin.de'
                    )
                ]
            ),
            'apmtools': NORTHTool(
                short_description='An example for analyzing atom probe data.',
                description='Miscellaneous tools from the atom probe community:\nCurrently APTyzer, paraprobe-toolbox, and APAV',
                image='gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-remote-tools-hub/apmtools-webtop',
                privileged=True,
                icon='jupyter_logo.svg',
                mount_path='/config',
                maintainer=[
                    NORTHToolMaintainer(
                        name='Markus Kühbach',
                        email='markus.kuehbach@physik.hu-berlin.de'
                    )
                ]
            ),
            'fiji': NORTHTool(
                short_description='ImageJ and Fiji for image processing',
                description='ImageJ and Fiji with amongst others several electron-microscopy specific plug-ins',
                image='gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-remote-tools-hub/fiji-webtop',
                privileged=True,
                icon='jupyter_logo.svg',
                mount_path='/config',
                maintainer=[
                    NORTHToolMaintainer(
                        name='Markus Kühbach',
                        email='markus.kuehbach@physik.hu-berlin.de'
                    )
                ]
            ),
            'frwr': NORTHTool(
                short_description='Inline electron holography by C. Koch',
                description='FRWR3 in-line holography/focus series reconstruction code',
                image='gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-remote-tools-hub/frwr-webtop',
                privileged=True,
                icon='jupyter_logo.svg',
                mount_path='/config',
                maintainer=[
                    NORTHToolMaintainer(
                        name='Markus Kühbach',
                        email='markus.kuehbach@physik.hu-berlin.de'
                    )
                ]
            ),
            'abtem': NORTHTool(
                short_description='Electronic structure supported image simulation for transmission electron microscopy.',
                description='VESTA, GPAW, and abTEM configured in one container for simulating images and diffraction patterns in transmission electron microscopy',
                image='gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-remote-tools-hub/abtem-webtop',
                privileged=True,
                icon='jupyter_logo.svg',
                mount_path='/config',
                maintainer=[
                    NORTHToolMaintainer(
                        name='Markus Kühbach',
                        email='markus.kuehbach@physik.hu-berlin.de'
                    )
                ]
            )
        }),
        description='The available north tools. Either the tools definitions as dict or a path to a .json file.')

    hub_service_api_token: str = Field('secret-token', description='''
        A secret token shared between NOMAD and the NORTH jupyterhub.
        This needs to be the token of an admin service.''')
