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

import os

import click
from tqdm import tqdm

from .cli import cli


@cli.command(
    help='Cleanse the given path by removing empty folders.',
    name='clean',
)
@click.option(
    '--path',
    type=str,
    help='Cleanse the given path by removing empty folders.',
)
def clean_staging(path):
    if 'staging' not in path:
        print('Path must contain "staging".')
        return

    print(f'Cleaning path: "{path}".')
    print('Are you sure you want to continue? (y/N)', end=' ')
    response = input()
    if response.lower() != 'y':
        print('Exiting...')
        return

    print('Cleaning...')

    def safe_remove(_p):
        try:
            os.rmdir(_p)
        except Exception:  # noqa
            pass

    for root, folders, _ in tqdm(os.walk(path, topdown=False)):
        for folder in folders:
            if not os.listdir(full_path := os.path.join(root, folder)):  # noqa
                safe_remove(full_path)
