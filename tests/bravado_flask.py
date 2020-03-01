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

from urllib.parse import urlencode
from bravado.http_client import HttpClient
from bravado.http_future import HttpFuture
from bravado_core.response import IncomingResponse
import json


class FlaskTestHttpClient(HttpClient):
    def __init__(self, flask_test_client, headers={}):
        self._flask_client = flask_test_client
        self._headers = headers

    def request(self, request_params, *args, **kwargs):
        '''
        Taken from `bravado.http_client.HttpClient`.

        Args:
            request_params (dict): complete request data. e.g. url, method, headers, body, params,
                connect_timeout, timeout, etc.
            operation (`bravado_core.operation.Operation`): operation that this http request
                is for. Defaults to None - in which case, we're obviously just retrieving a Swagger
                Spec.
            response_callbacks: List of callables to post-process the incoming response.
                Expects args incoming_response and operation.
            also_return_response: Consult the constructor documentation for
                `bravado.http_future.HttpFuture`.
        Returns:
            `bravado_core.http_future.HttpFuture`: HTTP Future object
        '''
        request_params.setdefault('headers', {}).update(self._headers)
        test_future = FlaskTestFutureAdapter(request_params, self._flask_client)

        return HttpFuture(test_future, FlaskTestResponseAdapter, *args, **kwargs)


class FlaskTestFutureAdapter:
    '''
    Mimics a :class:`concurrent.futures.Future` for the purposes of making it work with
    Bravado's :class:`bravado.http_future.HttpFuture` when simulating calls to a Falcon API.
    Those calls will be validated by Bravado.

    Args:
        request_params (dict): Request parameters provided to
            :class:`bravado.http_client.HttpClient` interface.
        falcon_api (`falcon.API`): API object to send the request to.
        response_encoding (str): Encoding that will be used to decode response's body.
            If set to None then the body won't be decoded.
    '''

    def __init__(self, request_params, flask_client, response_encoding='utf-8'):
        self._flask_client = flask_client
        self._request_params = request_params
        self._response_encoding = response_encoding

        self.timeout_errors = None
        self.connection_errors = None

    def result(self, **_):
        '''
        Args:
            **_: Ignore all the keyword arguments (right now it's just timeout) passed by Bravado.
        '''
        # Bravado will create the URL by appending request path to 'http://localhost'
        path = self._request_params['url'].replace('http://localhost', '')
        method = self._request_params.get('method')

        query = urlencode(self._request_params.get('params', {}), doseq=True)
        if query is not None and query != '':
            url = '%s?%s' % (path, query)
        else:
            url = path

        data = self._request_params.get('data')

        function = getattr(self._flask_client, method.lower())

        files = self._request_params.get('files', [])
        if len(files) > 1:
            raise NotImplementedError
        if len(files) == 1:
            _, (_, f) = files[0]
            data = f

        return function(
            url, headers=self._request_params.get('headers'), data=data)


class FlaskTestResponseAdapter(IncomingResponse):
    '''
    Wraps a response from Falcon test client to provide a uniform interface
    expected by Bravado's :class:`bravado.http_future.HttpFuture`.
    Args:
        flask_response: Response to a call simulated with flask's test client.
    '''

    def __init__(self, flask_response):
        self._response = flask_response

    @property
    def status_code(self):
        '''
        Returns:
            int: HTTP status code
        '''
        return self._response.status_code

    @property
    def text(self):
        '''
        Returns:
            str: Textual representation of the response's body.
        '''
        return self._response.data

    @property
    def reason(self):
        '''
        Returns:
            str: Reason-phrase of the HTTP response (e.g. "OK", or "Not Found")
        '''
        # status codes from Falcon look like this: "200 OK"
        return self._response.status[4:]

    @property
    def headers(self):
        '''
        Returns:
            dict: Headers attached to the response.
        '''
        return self._response.headers

    def json(self, **kwargs):
        '''
        Args:
            **kwargs: This is a part of the interface, but we don't do anything with it.
        Returns:
            dict: JSON representation of the response's body.
        '''
        return json.loads(self._response.data)
