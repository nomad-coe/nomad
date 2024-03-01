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

# This requires licensheaders. There is a bug in the released version and we need a
# manually pip install from a clone of git@github.com:johann-petrak/licenseheaders.git

import subprocess
import re
import os
import requests
import json
import yaml
import sys
import time


def str_presenter(dumper, data):
    if '\n' in data:  # check for multiline string
        return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')
    return dumper.represent_scalar('tag:yaml.org,2002:str', data)


yaml.add_representer(str, str_presenter)


def sh(cmd):
    result = subprocess.run(['bash', '-c', cmd], stdout=subprocess.PIPE)
    output = result.stdout.decode('utf-8')
    return output


def collect_submodules(filter: str = None):
    lines = sh('git submodule status').split('\n')
    for line in lines:
        line = line.strip()
        if line == '':
            continue

        commit, path, _ = line.split(' ', 2)

        if filter is not None:
            if re.search(filter, path):
                yield path, commit
        else:
            yield path, commit


def move_submodule_to_github(path: str):
    name = os.path.basename(path)
    github_url = f'https://github.com/nomad-coe/nomad-parser-{name}.git'
    github_giturl = f'git@github.com:nomad-coe/nomad-parser-{name}.git'

    if requests.get(github_url).status_code != 404:
        print(f'Skip {path}. It is already on GitHUB at {github_url}.')
        return

    origin = next(
        line for line in sh('git remote -v').split('\n') if line.startswith('origin')
    )

    origin = re.split(r'\t| ', origin)[1]
    if 'github.com/nomad-coe' in origin:
        # print(f'Skip {path}. It already has GitHUB origin {origin}.')
        return

    print(f'Moving {path} now.')

    # get full history from current remove
    sh('git fetch')
    sh('git fetch --tags')

    # checkout master branch
    sh('git checkout origin/master')

    # make changes
    # 1. license file
    license_file = os.path.join(os.path.dirname(__file__), '../../../LICENSE')
    sh(f'cp {license_file} .')
    if os.path.exists('LICENSE.txt'):
        os.remove('LICENSE.txt')
    # 2. update copyright headers
    header_templ = os.path.join(os.path.dirname(__file__), 'apache-2.tmpl')
    sh(
        f'licenseheaders -t {header_templ} -o "The NOMAD Authors" -u "https://nomad-lab.eu" -n NOMAD -x "*.yaml" "*.yml"'
    )
    # 3. create an authors file
    authors_file = 'AUTHORS'
    authors = [
        re.split(r'\t| ', line.strip(), 1)[1]
        for line in sh('git shortlog -sne').splitlines()
    ]

    with open(authors_file, 'wt') as f:
        f.write('\n'.join(authors))
        f.write('\n')
    # 4. update the metadata
    with open('metadata.yaml', 'rt') as f:
        parser_metadata = yaml.load(f, Loader=yaml.SafeLoader)
        title = parser_metadata['codeLabel']
        parser_metadata['parserGitUrl'] = github_url
        parser_metadata['preamble'] = ''
        parser_metadata['codeName'] = name
    with open('metadata.yaml', 'wt') as f:
        yaml.dump(parser_metadata, f)
    sh(f'nomad dev update-parser-readmes --parser {name}')

    # commit changes
    sh('git add -A')
    sh('git commit -a -m "Prepared license and copyright headers."')

    # create a git-repo at GitHUB
    result = requests.post(
        'https://api.github.com/orgs/nomad-coe/repos',
        headers={
            'Accept': 'application/vnd.github.v3+json',
            'Authorization': f'token {os.environ["GITHUB_TOKEN"]}',
        },
        data=json.dumps(
            {
                'name': f'nomad-parser-{name}',
                'description': f"This is a NOMAD parser for {title}. It will read {title} input and output files and provide all information in NOMAD's unified Metainfo based Archive format.",
                'private': False,
                'visibility': True,
                'has_issues': True,
                'has_projects': False,
                'has_wiki': False,
                'team_id': 4443059,
            }
        ),
    )

    if result.status_code >= 400:
        print('Could not create github repository: ', result.text)
        return

    time.sleep(3)

    # set team permission
    result = requests.put(
        f'https://api.github.com/orgs/nomad-coe/teams/nomad-parser-developer/repos/nomad-coe/nomad-parser-{name}',
        headers={
            'Accept': 'application/vnd.github.v3+json',
            'Authorization': f'token {os.environ["GITHUB_TOKEN"]}',
        },
        data=json.dumps({'permission': 'push'}),
    )

    if result.status_code >= 400:
        print('Could not set tream permission repository: ', result.text)
        return

    # reset remote
    sh(f'git remote set-url origin {github_giturl}')

    # push
    sh(f'git branch -D master')
    sh(f'git checkout -b master')
    sh(f'git push --all')
    sh(f'git push --tags')

    # change submodule
    nomad_dir = os.path.join(os.path.dirname(__file__), '../../..')
    os.chdir(nomad_dir)
    sh(f'git config --file=.gitmodules submodule.{path}.url {github_url}')
    sh(f'git config --file=.gitmodules submodule.{path}.branch master')


def run_move_submodules_to_github():
    working_dir = os.path.abspath(os.curdir)

    parser_submodules = sys.argv[1:]

    if len(parser_submodules) == 0:
        parser_submodules = [path for path, _ in collect_submodules(filter=r'parsers')]

    print(parser_submodules)
    for submodule in parser_submodules:
        path = submodule
        os.chdir(path)
        move_submodule_to_github(path)
        os.chdir(working_dir)


if __name__ == '__main__':
    global __file__
    __file__ = os.path.abspath(__file__)
    run_move_submodules_to_github()
