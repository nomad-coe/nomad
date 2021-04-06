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

from typing import Dict, Iterator, Any
from types import FunctionType
import urllib
import sys
import inspect
from fastapi import Request, Query, HTTPException  # pylint: disable=unused-import
from pydantic import ValidationError, BaseModel  # pylint: disable=unused-import
import zipstream

if sys.version_info >= (3, 7):
    import zipfile
else:
    import zipfile37 as zipfile  # pragma: no cover


def parameter_dependency_from_model(name: str, model_cls):
    '''
    Takes a pydantic model class as input and creates a dependency with corresponding
    Query parameter definitions that can be used for GET
    requests.

    This will only work, if the fields defined in the input model can be turned into
    suitable query parameters. Otherwise fastapi will complain down the road.

    Arguments:
        name: Name for the dependency function.
        model_cls: A ``BaseModel`` inheriting model class as input.
    '''
    names = []
    annotations: Dict[str, type] = {}
    defaults = []
    for field_model in model_cls.__fields__.values():
        field_info = field_model.field_info

        names.append(field_model.name)
        annotations[field_model.name] = field_model.outer_type_
        defaults.append(Query(field_model.default, description=field_info.description))

    code = inspect.cleandoc('''
    def %s(%s):
        try:
            return %s(%s)
        except ValidationError as e:
            errors = e.errors()
            for error in errors:
                error['loc'] = ['query'] + list(error['loc'])
            raise HTTPException(422, detail=errors)

    ''' % (
        name, ', '.join(names), model_cls.__name__,
        ', '.join(['%s=%s' % (name, name) for name in names])))

    compiled = compile(code, 'string', 'exec')
    env = {model_cls.__name__: model_cls}
    env.update(**globals())
    func = FunctionType(compiled.co_consts[0], env, name)
    func.__annotations__ = annotations
    func.__defaults__ = (*defaults,)

    return func


class File(BaseModel):
    path: str
    f: Any
    size: int


def create_streamed_zipfile(
        files: Iterator[File],
        compress: bool = False) -> Iterator[bytes]:

    '''
    Creates a streaming zipfile object that can be used in fastapi's ``StreamingResponse``.
    '''

    def path_to_write_generator():
        for file_obj in files:
            def content_generator():
                while True:
                    data = file_obj.f.read(1024 * 64)
                    if not data:
                        break
                    yield data

            yield dict(
                arcname=file_obj.path,
                iterable=content_generator(),
                buffer_size=file_obj.size)

    compression = zipfile.ZIP_DEFLATED if compress else zipfile.ZIP_STORED
    zip_stream = zipstream.ZipFile(mode='w', compression=compression, allowZip64=True)
    zip_stream.paths_to_write = path_to_write_generator()

    for chunk in zip_stream:
        yield chunk


def create_responses(*args):
    return {
        status_code: response
        for status_code, response in args}


def update_url_query_arguments(original_url: str, **kwargs) -> str:
    '''
    Takes an url, and returns a new url, obtained by updating the query arguments in the
    `original_url` as specified by the kwargs.
    '''
    scheme, netloc, path, params, query, fragment = urllib.parse.urlparse(original_url)
    query_dict = urllib.parse.parse_qs(query)
    for k, v in kwargs.items():
        if v is None:
            # Unset the value
            if k in query_dict:
                query_dict.pop(k)
        else:
            query_dict[k] = [str(v)]
    query = urllib.parse.urlencode(query_dict, doseq=True)
    return urllib.parse.urlunparse((scheme, netloc, path, params, query, fragment))
