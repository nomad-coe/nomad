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

from celery import Celery, group, subtask
import re
import logging
import nomad.config as config
import nomad.files as files

broker_url = 'pyamqp://%s:%s@localhost//' % (config.rabbitmq.user, config.rabbitmq.password)
backend_url = 'rpc://localhost'
app = Celery('nomad.processing', backend=backend_url, broker=broker_url)
app.conf.update(
    accept_content=['pickle'],
    task_serializer='pickle',
    result_serializer='pickle',
)

LOGGER = logging.getLogger(__name__)


class Parser():
    """
    Instances specify a parser. It allows to find *main files* from  given uploaded
    and extracted files. Further, allows to run the parser on those 'main files'.
    """
    def __init__(self, name, main_file_re, main_contents_re):
        self.name = name
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

    def run(self, upload, filename):
        pass


parsers = [
    Parser(
        name='VaspRun',
        main_file_re=r'^.*\.xml$',
        main_contents_re=(
            r'^\s*<\?xml version="1\.0" encoding="ISO-8859-1"\?>\s*'
            r'?\s*<modeling>'
            r'?\s*<generator>'
            r'?\s*<i name="program" type="string">\s*vasp\s*</i>'
            r'?'
        )
    )
]

parser_dict = {parser.name: parser for parser in parsers}


@app.task()
def find_mainfiles(upload):
    mainfile_specs = list()
    for filename in upload.filelist:
        for parser in parsers:
            if parser.is_mainfile(upload, filename):
                mainfile_specs.append((upload, filename, parser.name))

    return mainfile_specs


@app.task()
def open_upload(upload_id):
    upload = files.upload(upload_id)
    upload.open()
    return upload


@app.task()
def close_upload(upload):
    upload.close()


@app.task()
def parse(mainfile_spec):
    upload, mainfile, parser = mainfile_spec
    LOGGER.debug('Start parsing mainfile %s/%s with %s.' % (upload, mainfile, parser))
    parser_dict[parser].run(upload, mainfile)

    return True


@app.task()
def dmap(it, callback):
    callback = subtask(callback)
    return group(callback.clone([arg, ]) for arg in it)()

if __name__ == '__main__':
    upload_id = 'examples_vasp.zip'
    parsing_workflow = (open_upload.s(upload_id) | find_mainfiles.s() | dmap.s(parse.s()))

    print(~parsing_workflow)
