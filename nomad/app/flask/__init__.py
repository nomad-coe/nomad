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

'''
This module comprises the nomad@FAIRDI APIs. Currently there is NOMAD's official api, optimade api,
and dcat api. The app module also servers documentation, gui, and
alive.
'''
from flask import Flask
from flask_cors import CORS

from .docs import blueprint as docs_blueprint
from .dist import blueprint as dist_blueprint
from .gui import blueprint as gui_blueprint


app = Flask(__name__)
''' The Flask app that serves all APIs. '''

CORS(app)

app.register_blueprint(docs_blueprint, url_prefix='/docs')
app.register_blueprint(dist_blueprint, url_prefix='/dist')
app.register_blueprint(gui_blueprint, url_prefix='/gui')


@app.route('/alive')
def alive():
    ''' Simple endpoint to utilize kubernetes liveness/readiness probing. '''
    return "I am, alive!"
