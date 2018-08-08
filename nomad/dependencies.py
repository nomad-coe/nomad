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
Integration of nomad projects into the processing
=================================================

Parsers are developed as independed, individual python programs in their own GIT repositories.
They are build on a common modules called *python-common*, also in a separate GIT.
All parsers depend on the *meta-info*, which is also maintained in its own GIT.

Assumption about parsers
------------------------
For now, we make a few assumption about parsers
- they always work on the same *meta-info* version
- they have no conflicting python requirments
- they can be loaded at the same time and can be used within the same python process
- they are uniquely identified by a GIT URL and publicly accessible
- their version is uniquly identified by a GIT commit SHA

Preparing dependencies and parsers
----------------------------------

To make GIT maintained python modules available, we use:

.. autoclass:: nomad.parsers.PythonGit

Parsers, as a special case for a GIT maintained python modules, can be used via:

.. autoclass:: nomad.parsers.Parser
"""
import re
import sys
import os
import os.path
from git import Repo, Git
try:
    from pip import main as pip
except:
    from pip._internal import main as pip
import importlib
import logging
import subprocess

_meta_info_path = './submodules/nomad-meta-info/meta_info/nomad_meta_info/'
_logger = logging.getLogger(__name__)
base_dir = './.dependencies'


class PythonGitError(Exception):
    def __init__(self, msg, repo):
        msg = '%s [%s]' % (msg, repo)
        super().__init__(msg)


class PythonGit():
    """Represents a python module in a git repository.

    It allows to fetch a specific commit, install all requirements to
    the current python environment, and check the installation via module import.

    This is only useful before you want to use the respective module in a different
    python process, because it will not try to reload any already loaded modules into
    the current python process.
    """
    def __init__(self, name, git_url, git_commit):
        """
        Args:
            name: A name that determines the download path, can contain '/' for sub dirs.
                  Names are important, because modules might use relatives paths between
                  them.
            git_url: A publically available and fetchable url to the GIT repository.
            git_commit: The full commit SHA of the desired commit.
        """
        super().__init__()
        self.name = name
        self.git_url = git_url
        self.git_commit = git_commit

    def _run_pip_install(self, *args):
        pipcode = 0

        # some weird interaction of pip and virtualenv causes a bug that does
        # not allow to install in docker due to a wrong PIP_REQ_TRACKER path. This
        # is a workarround.
        pip_req_tracker_key = 'PIP_REQ_TRACKER'
        env = dict(os.environ)
        if pip_req_tracker_key in env:
            del(env['PIP_REQ_TRACKER'])
        pipcode = subprocess.call(
            [sys.executable, '-m', 'pip', 'install'] + list(args),
            env=env)

        if pipcode != 0:
            raise PythonGitError(
                'Could not install (pip return code=%s)' % pipcode, repo=self)

    def prepare(self):
        """Makes sure that the repository is fetched, at the right commit, and installed.

        Raises:
            PythonGitError: if something went wrong.
        """
        # check/change working directory
        old_cwd = os.getcwd()
        try:
            cwd = os.path.join(base_dir, self.name)
            if not os.path.exists(cwd):
                os.makedirs(cwd)
            os.chdir(cwd)

            _logger.info('check git/do init with origin %s for %s' % (self.git_url, self.name))
            if os.path.exists('.git'):
                git = Repo('./')
            else:
                git_cmd = Git('./')
                git_cmd.init()
                git = Repo('./')
                origin = git.create_remote('origin', self.git_url)

            _logger.info('pull %s for %s' % (self.git_commit, self.name))
            origin = git.remote('origin')
            origin.pull(self.git_commit)

            if os.path.exists('requirements.txt'):
                _logger.info('install requirements.txt for %s' % self.name)
                self._run_pip_install('-r', 'requirements.txt')

            if os.path.exists('setup.py'):
                _logger.info('install setup.py for %s' % self.name)
                self._run_pip_install('-e', '.')

        except PythonGitError as e:
            raise e
        except Exception as e:
            raise PythonGitError(
                'Unexpected exception during preparation: %s' % e, repo=self)
        finally:
            os.chdir(old_cwd)
        pass

    def __repr__(self):
        return self.name


class Parser():
    """
    Instances specify a parser. It allows to find *main files* from  given uploaded
    and extracted files. Further, allows to run the parser on those 'main files'.
    """
    def __init__(self, python_git, parser_class_name, main_file_re, main_contents_re):
        self.name = python_git.name
        self.python_git = python_git
        self.parser_class_name = parser_class_name
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
        from nomadcore.parser_backend import JsonParseEventsWriterBackend
        module_name = self.parser_class_name.split('.')[:-1]
        parser_class = self.parser_class_name.split('.')[1]
        module = importlib.import_module('.'.join(module_name))
        Parser = getattr(module, parser_class)
        parser = Parser(backend=JsonParseEventsWriterBackend)
        parser.parse(mainfile)


class VASPRunParser(Parser):
    def __init__(self):
        super().__init__(
            python_git=PythonGit(
                name='parsers/vasp',
                git_url='https://gitlab.mpcdf.mpg.de/nomad-lab/parser-vasp.git',
                git_commit='nomad-xt'),
            parser_class_name='vaspparser.VASPParser',
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

others = [
    PythonGit(
        name='nomad-meta-info',
        git_url='https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info.git',
        git_commit='1.6.0'),
    PythonGit(
        name='python_common',
        git_url='https://gitlab.mpcdf.mpg.de/nomad-lab/python-common.git',
        git_commit='master'),  # COMMIT b6f64de2149f95da6a79c4f86fd909b1dcfc23e8
]


def prepare():
    for python_git in others:
        python_git.prepare()

    for parser in parsers:
        parser.python_git.prepare()

if __name__ == '__main__':
    _logger.setLevel(logging.DEBUG)
    prepare()
