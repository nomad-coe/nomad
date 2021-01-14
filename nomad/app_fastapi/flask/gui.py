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

from flask import Blueprint, request
import os.path

gui_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), 'static/gui'))
configuired_gui_folder = os.path.join(gui_folder, '../.gui_configured')
if os.path.exists(configuired_gui_folder):
    gui_folder = configuired_gui_folder
blueprint = Blueprint('gui', __name__, static_url_path='/', static_folder=gui_folder)


@blueprint.after_request
def add_header(response):
    if request.url.endswith('index.html'):
        response.headers['Cache-Control'] = 'no-cache, no-store, must-revalidate'
    if request.url.endswith('.js'):
        response.headers['Cache-Control'] = 'no-cache, must-revalidate'
    return response
