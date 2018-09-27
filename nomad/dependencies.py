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
This module allows to configure and install all necessary legecy nomad GIT repositories
to process (parser, normalizer, etc.) uploaded calculations.

Parsers are developed as independed, individual python programs in their own GIT repositories.
They are build on a common modules called *python-common*, also in a separate GIT.
All parsers depend on the *meta-info*, which is also maintained in its own GIT.

Preparing dependencies
----------------------

To make GIT maintained python modules available, we use:

.. autoclass:: PythonGit
    :members:


Dependencies are configured in

.. autodata:: dependencies


To install all dependencies use

.. autofunction:: prepare
"""
import sys
import os
import os.path
import logging
import subprocess
import shutil

_meta_info_path = './submodules/nomad-meta-info/meta_info/nomad_meta_info/'
_logger = logging.getLogger(__name__)
base_dir = './.dependencies'


class PythonGitError(Exception):
    def __init__(self, msg, repo):
        msg = '%s [%s]' % (msg, repo)
        super().__init__(msg)


class PythonGit():
    """
    Represents a python module in a git repository. It allows to fetch a specific commit,
    install all requirements to the current python environment, and check the installation
    via module import.

    This is only useful before you want to use the respective module in a different
    python process, because it will not try to reload any already loaded modules into
    the current python process.

    Arguments:
        name: A name that determines the download path, can contain '/' for sub dirs.
                Names are important, because modules might use relatives paths between
                them.
        git_url: A publically available and fetchable url to the GIT repository.
        git_commit: The full commit SHA of the desired commit.
    """
    def __init__(self, name: str, git_url: str, git_commit: str) -> None:
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

    def prepare(self, dev: bool = False) -> None:
        """
        Makes sure that the repository is fetched, at the right commit, and installed.

        Arguments:
            dev (bool): Indicate dev install (uses pip with -e). Default is False.

        Raises:
            PythonGitError: if something went wrong.
        """
        # import late because git will not be available in production
        from git import Repo, Git

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
            origin.pull(self.git_commit, depth=1)

            if os.path.exists('requirements.txt'):
                _logger.info('install requirements.txt for %s' % self.name)
                self._run_pip_install('-r', 'requirements.txt')

            if os.path.exists('setup.py'):
                _logger.info('install setup.py for %s' % self.name)
                if dev:
                    self._run_pip_install('-e', '.')
                else:
                    self._run_pip_install('.')

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


dependencies = [
    PythonGit(
        name='nomad-meta-info',
        git_url='https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info.git',
        git_commit='nomad-fair'),
    PythonGit(
        name='python_common',
        git_url='https://gitlab.mpcdf.mpg.de/nomad-lab/python-common.git',
        git_commit='nomad-fair'),
    PythonGit(
        name='parsers/vasp',
        git_url='https://gitlab.mpcdf.mpg.de/nomad-lab/parser-vasp.git',
        git_commit='nomad-fair'),
    PythonGit(
        name='parsers/exciting',
        git_url='https://gitlab.mpcdf.mpg.de/nomad-lab/parser-exciting.git',
        git_commit='nomad-fair'),
    PythonGit(
        name='parsers/fhi-aims',
        git_url='https://gitlab.mpcdf.mpg.de/nomad-lab/parser-fhi-aims.git',
        git_commit='nomad-fair'),
    PythonGit(
        name='normalizers/stats',
        git_url='https://gitlab.mpcdf.mpg.de/nomad-lab/normalizer-stats.git',
        git_commit='nomad-fair'),
    PythonGit(
        name='normalizers/symmetry',
        git_url='https://gitlab.mpcdf.mpg.de/nomad-lab/normalizer-symmetry',
        git_commit='nomad-fair'),
    PythonGit(
        name='normalizers/system-type',
        git_url='https://gitlab.mpcdf.mpg.de/nomad-lab/normalizer-system-type',
        git_commit='nomad-fair')
]

dependencies_dict = {dependency.name: dependency for dependency in dependencies}


def prepare(*args, **kwargs) -> None:
    """
    Installs all dependencies from :data:`dependencies` and :data:`parsers`.
    """
    for python_git in dependencies:
        python_git.prepare(*args, **kwargs)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Install dependencies from NOMAD-coe.')
    parser.add_argument('--dev', help='pip install with -e', action='store_true')
    parser.add_argument('--clean', help='remove old deps first', action='store_true')

    args = parser.parse_args()

    _logger.setLevel(logging.DEBUG)
    if args.clean and os.path.exists(base_dir):
        shutil.rmtree(base_dir)

    prepare(dev=args.dev)
