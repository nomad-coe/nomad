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

import click
import json

from nomad import infrastructure, datamodel

from .admin import admin


@admin.command(help='Import users to keycloak from a JSON file.', name='import')
@click.argument('PATH_TO_USERS_FILE', type=str, nargs=1)
def import_command(path_to_users_file):
    with open(path_to_users_file, 'rt') as f:
        users = json.load(f)

    infrastructure.setup_logging()

    for user_dict in users:
        password = user_dict.pop('password')
        user = datamodel.User(**user_dict)
        infrastructure.keycloak.add_user(user, bcrypt_password=password, invite=False)
        print('Imported %s' % user.name)
