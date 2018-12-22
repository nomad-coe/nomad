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

from werkzeug.wsgi import DispatcherMiddleware

from nomad.api import app
from nomad import config


def run_dev_server(*args, **kwargs):
    def simple(env, resp):
        resp(b'200 OK', [(b'Content-Type', b'text/plain')])
        return [
            ('Development nomad api server. Api is served under %s/.' %
                config.services.api_base_path).encode('utf-8')]

    app.wsgi_app = DispatcherMiddleware(simple, {config.services.api_base_path: app.wsgi_app})
    app.run(*args, **kwargs)


if __name__ == '__main__':
    run_dev_server(debug=True, port=8000)
