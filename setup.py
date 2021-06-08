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

from setuptools import setup, find_packages
from subprocess import call
from setuptools.command.install import install as setup_install
import runpy
import os
import sys
import json
import re
import fastentrypoints  # pylint: disable=unused-import

'''
This setup.py works differently for creating a distribution than installing for
development. The idea is that for distributions, we compile the requirements, packages, and
data files from all dependencies and install everything under one package called 'nomad-lab'.
For development we install each dependency as its own project. This is mainly due to
pip's limitation of not being able to install from multiple source folders in develop
mode.

The compilation process can by run with ``python setup.py compile``. It will create a
file with setup kwargs calles ``./setup.json``. If this file exists setup will be called
with compiled args. If this file does not exist, setup will be called with regular args.

This setup makes use of ``extras_require``. The extras are parsing, infrastructure, dev,
and all. Where all is the union of the other ones. The extras are parsed from the
requirements.txt where specific comments are used to assign an extra to requirements.
'''


ignore_extra_requires = ['optimade']
''' Dependencies where the extra_requires should not be added '''


def parse_requirements():
    '''
    Parses the requirements.txt file to extras install and extra requirements.
    Sections headed with # [extra] are assigned to the 'extra' extra.

    Returns:
        Tuple with install and extra requires passible to :func:`setuptools.setup`.
    '''
    with open('requirements.txt', 'rt') as f:
        lines = f.readlines()

    extras_require = {}
    requires = []
    all_requires = []
    current = None
    for line in lines:
        line = line.strip()

        if line == '':
            continue

        match = re.match(r'^#\s*\[([a-zA-Z0-9_]+)\]$', line)
        if match:
            extra = match.group(1)
            current = list()
            extras_require[extra] = current
        elif line.startswith('#'):
            continue
        else:
            line = line.split('#')[0].strip()
            if current is None:
                requires.append(line)
            else:
                current.append(line)
                all_requires.append(line)

    extras_require['all'] = all_requires

    return requires, extras_require


class install(setup_install):
    def __post_install(self, dir):
        if os.name == 'posix':
            call(['./auto_complete_install.sh'])

    def run(self):
        setup_install.run(self)
        self.execute(self.__post_install, (self.install_lib, ), msg='installing autocompletion')


def compile_dependency_setup_kwargs(paths, **kwargs):
    import setuptools
    import distutils.core as distutils_core

    # collect all kwargs from all setup.pys
    results = {}
    current = {}

    def patched_setup(*args, **kwargs):
        assert len(args) == 0
        assert current['name'] not in results, 'current is %s' % current['name']

        results[current['name']] = {
            'meta': dict(**current),
            'kwargs': kwargs
        }

    if len(kwargs) > 0:
        current['name'] = kwargs.get('name', 'nomad')
        current['directory'] = './'
        current['setup.py'] = './setup.py'
        patched_setup(**kwargs)

    setuptoolss = [setuptools, distutils_core]
    orig_setups = []
    for st in setuptoolss:
        orig_setups.append((st, getattr(st, 'setup')))
        setattr(st, 'setup', patched_setup)

    for path in paths:
        for root, _, files in os.walk(path):
            for file in files:
                if os.path.basename(file) == 'setup.py':
                    setup_path = os.path.join(root, file)
                    current['name'] = os.path.basename(os.path.dirname(setup_path))
                    current['directory'] = os.path.dirname(setup_path)
                    current['setup.py'] = setup_path
                    cwd = os.getcwd()
                    os.chdir(os.path.dirname(setup_path))
                    try:
                        runpy.run_path(file, run_name='__main__')
                    except Exception:
                        import traceback
                        traceback.print_exc()
                        print('Could not run %s' % setup_path)
                        sys.exit(1)
                    finally:
                        os.chdir(cwd)

    for st, setup in orig_setups:
        setattr(st, 'setup', setup)

    # combine the kwargs
    all_packages = []
    all_package_dir = {}
    all_package_data = {}
    all_install_requires = {}
    all_names = set()
    for _, setup_data in results.items():
        meta = setup_data['meta']
        local_kwargs = setup_data['kwargs']

        if 'name' in local_kwargs:
            all_names.add(local_kwargs['name'])

        name = local_kwargs.get('name', meta['name'])

        # 1. packages, package_dir
        package_dir = local_kwargs.get('package_dir', {'': './'})
        packages = local_kwargs.get('packages', [])

        assert len(package_dir) == 1
        assert '' in package_dir

        for package in packages:
            root = package.split('.')[0]
            if root not in package_dir and root not in all_package_dir:
                all_package_dir[root] = os.path.normpath(
                    os.path.join(meta['directory'], package_dir[''], root))
            all_packages.append(package)

        # 2. package_data
        all_package_data.update(**local_kwargs.get('package_data', {}))

        # 3. requires
        local_install_requires = set()
        if name not in ignore_extra_requires:
            for extra_require in local_kwargs.get('extras_require', {}).values():
                for require in extra_require:
                    local_install_requires.add(require)

        for require in local_kwargs.get('install_requires', []):
            local_install_requires.add(require)
        all_install_requires[name] = local_install_requires

    # automatically add parser deps
    for _, setup_data in results.items():
        if 'parsers' in setup_data['meta']['setup.py']:
            parsing = kwargs.setdefault('extras_require', {}).setdefault('parsing', [])
            all = kwargs.setdefault('extras_require', {}).setdefault('all', [])
            for require in setup_data['kwargs'].get('install_requires', []):
                if require not in parsing:
                    parsing.append(require)
                    all.append(require)

    def replace_own_packages(requires):
        ''' replaces nomad dependencies with their requirements '''
        for other in all_names:
            if other in requires:
                requires.remove(other)
                requires.extend(all_install_requires[other])

        # remove dups
        sorted_requires = sorted(requires)
        sorted_normalized = [re.match(r'^([a-zA-Z0-9_\.\-]*).*$', r).group(1) for r in sorted_requires]
        to_remove = []
        for i, require in enumerate(sorted_normalized):
            if i + 1 < len(sorted_requires) and require == sorted_normalized[i + 1]:
                to_remove.append(i)
        for i in to_remove:
            requires.remove(sorted_requires[i])

    # run twice because dependencies can be dependencies of dependencies
    for _ in range(0, 2):
        for extra_require in kwargs.get('extras_require', {}).values():
            replace_own_packages(extra_require)
        replace_own_packages(kwargs.get('install_requires', []))

    kwargs.update(**{
        'package_dir': all_package_dir,
        'packages': all_packages,
        'package_data': all_package_data
    })

    return kwargs


def setup_kwargs():
    from nomad import config

    install_requires, extras_require = parse_requirements()
    with open("README.md", "r") as fh:
        long_description = fh.read()

    return dict(
        name='nomad-lab',
        author='NOMAD Laboratory',
        author_email='markus.scheidgen@physik.hu-berlin.de',
        url='https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR',
        version=config.meta.version,
        license='APACHE 2.0',
        description='The NOvel MAterials Discovery (NOMAD) Python package',
        long_description=long_description,
        long_description_content_type="text/markdown",
        package_dir={'': './'},
        packages=['nomad.%s' % pkg for pkg in find_packages('./nomad')] + ['nomad'],
        setup_requires=['pip', 'setuptools', 'wheel', 'fastentrypoints', 'numpy', 'pyyaml'],
        install_requires=install_requires,
        extras_require=extras_require,
        include_package_data=True,
        python_requires='>=3.6',
        entry_points='''
            [console_scripts]
            nomad=nomad.cli:run_cli
        ''')


if __name__ == '__main__':
    if len(sys.argv) == 2 and sys.argv[1] == 'compile':
        kwargs = compile_dependency_setup_kwargs(['dependencies'], **setup_kwargs())
        kwargs['package_data']['optimade.grammar'] = ['*.lark']
        with open('setup.json', 'wt') as f:
            json.dump(kwargs, f, indent=2)
        sys.exit(0)

    if len(sys.argv) == 2 and sys.argv[1] == 'dry':
        import pprint
        pprint.pprint(
            compile_dependency_setup_kwargs(['dependencies'], **setup_kwargs()))
        sys.exit(0)

    if os.path.exists('setup.json'):
        with open('setup.json', 'rt') as f:
            kwargs = json.load(f)

    else:
        kwargs = setup_kwargs()

    setup(cmdclass={'install': install}, **kwargs)
