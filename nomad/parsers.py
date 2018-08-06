# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Integration of parsers into the processing
==========================================

Parsers are developed as independed, individual python programs in their own GIT repositories.
They are build on a common modules called *python-common*, also in a separate GIT.
All parsers depend on the *meta-info*, which is also maintained in its own GIT.

Assumption about parsers
------------------------
For now, we make a few assumption about parsers
- they always work on the same *meta-inf*
- they have no conflicting python requirments
- they can be loaded at the same time and can be used within the same python process
- they are uniquely identified by a GIT URL and publicly accessible
- their version is uniquly identified by a GIT commit SHA

Preparing dependencies and parsers during python run-time
---------------------------------------------------------
To make GIT maintained python modules available, we use:

.. autoclass:: nomad.parsers.PythonGitRepository

Parsers, as a special case for a GIT maintained python modules, can be used via:

.. autoclass:: nomad.parsers.Parser
"""
import re
import os
import os.path
from git import Repo, Git
try:
    from pip import main as pip
except:
    from pip._internal import main as pip
import importlib

from nomadcore.parser_backend import JsonParseEventsWriterBackend

_meta_info_path = './submodules/nomad-meta-info/meta_info/nomad_meta_info/'

base_dir = './.dependencies'


class PythonGitRepositoryError(Exception):
    def __init__(self, msg, repo):
        msg = '%s [%s]' % (msg, repo)
        super().__init__(msg)


class PythonGitRepository():
    """Represents a python module in a git repository.

    It allows to fetch a specific commit, install all requirements to
    the current python environment, and check the installation via module import.
    """
    def __init__(self, name, git_url, git_commit, modules=[]):
        """
        Args:
            name: A name that determines the download path, can contain '/' for sub dirs.
            git_url: A publically available and fetchable url to the GIT repository.
            git_commit: The full commit SHA of the desired commit.
            modules: A list of python module names that is used to confirm the installation.
        """
        super().__init__()
        self.name = name
        self.git_url = git_url
        self.git_commit = git_commit
        self.modules = modules

    def prepare(self, force_install=False):
        """Makes sure that the repository is fetched, at the right commit, and installed.

        Args:
            force_install: default is *False*. Allows to force install, e.g. after git commit or
                url change.

        Raises:
            PythonGitRepositoryError: if something went wrong.
        """
        # check/change working directory
        old_cwd = os.getcwd()
        try:
            cwd = os.path.join(base_dir, self.name)
            if not os.path.exists(cwd):
                os.makedirs(cwd)
            os.chdir(cwd)

            # check git/do init
            if os.path.exists('.git'):
                git = Repo('./')
            else:
                git_cmd = Git('./')
                git_cmd.init()
                git = Repo('./')
                origin = git.create_remote('origin', self.git_url)

            # check commit/checkout
            if 'master' not in git.heads:
                origin = git.remote('origin')
                origin.fetch(self.git_commit)
                git.create_head('master', self.git_commit)
            elif self.git_commit != git.heads.master.commit:
                origin = git.remote('origin')
                origin.fetch(self.git_commit)
            assert self.git_commit != git.heads.master.commit, \
                'Actual and desired commit do not match'
            git.heads.master.checkout()

            # check install
            def is_installed():
                for module in self.modules:
                    module_spec = importlib.util.find_spec(module)
                    if module_spec is None:
                        return False
                return True
            if is_installed() and not force_install:
                return

            # check/install requirements.txt
            if os.path.exists('requirements.txt'):
                # try twice to support circular dependencies
                for _ in range(1, 2):
                    pipcode = pip(['install', '-r', 'requirements.txt'])
                    if pipcode == 0:
                        break
                if pipcode != 0:
                    raise PythonGitRepositoryError(
                        'Could not install requirements (pip code=%s)' % pipcode, self)

            # check/install setup.py
            if os.path.exists('setup.py'):
                pipcode = pip(['install', '-e', '.'])
                if pipcode != 0:
                    raise PythonGitRepositoryError(
                        'Could not install (pip code=%s)' % pipcode, repo=self)

            # check install again
            if not is_installed():
                raise PythonGitRepositoryError(
                    'Some modules are not installed after install', repo=self)

            # reload, loaded modules when installed because of force_install
            # TODO
        except PythonGitRepositoryError as e:
            raise e
        except Exception as e:
            raise PythonGitRepositoryError(
                'Unexpected exception during preparation: %s' % e, repo=self)
        finally:
            os.chdir(old_cwd)
        pass


class Parser(PythonGitRepository):
    """
    Instances specify a parser. It allows to find *main files* from  given uploaded
    and extracted files. Further, allows to run the parser on those 'main files'.
    """
    def __init__(self, name, git_url, git_commit, parser, main_file_re, main_contents_re):
        modules = ['.'.join(parser.split('.')[:-1])]
        super().__init__(
            os.path.join('parsers', name), git_url, git_commit, modules=modules)
        self.parser = parser
        self._main_file_re = re.compile(main_file_re)
        self._main_contents_re = re.compile(main_contents_re)

    def is_mainfile(self, upload, filename):
        if self._main_file_re.match(filename):
            file = None
            try:
                file = upload.open_file(filename)
                return self._main_contents_re.match(file.read(500))
            finally:
                if file:
                    file.close()

    def run(self, mainfile):
        module_name = self.parser.split('.')[:-1]
        parser_class = self.parser.split('.')[1]
        module = importlib.import_module('.'.join(module_name))
        Parser = getattr(module, parser_class)
        parser = Parser(backend=JsonParseEventsWriterBackend)
        parser.parse(mainfile)


class VASPRunParser(Parser):
    def __init__(self):
        super().__init__(
            name='VASPRunParser',
            git_url='git@gitlab.mpcdf.mpg.de:nomad-lab/parser-vasp.git',
            git_commit='ddf8495944fbbcb62801f69b2c2c6c3d6099129d',
            parser='vaspparser.VASPParser',
            main_file_re=r'^.*\.xml$',
            main_contents_re=(
                r'^\s*<\?xml version="1\.0" encoding="ISO-8859-1"\?>\s*'
                r'?\s*<modeling>'
                r'?\s*<generator>'
                r'?\s*<i name="program" type="string">\s*vasp\s*</i>'
                r'?')
        )

parsers = [
    VASPRunParser()
]
parser_dict = {parser.name: parser for parser in parsers}


def prepare_parsers(force_install=False):
    for parser in parsers:
        parser.prepare(force_install=force_install)


if __name__ == '__main__':
    prepare_parsers(force_install=True)
