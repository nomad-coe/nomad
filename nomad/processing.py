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

from celery import Celery, chain, chord, group
import re
import nomad.config as config
import nomad.files as files

broker_url = 'pyamqp://%s:%s@localhost//' % (config.rabbitmq.user, config.rabbitmq.password)
backend_url = 'rpc://localhost'
app = Celery('nomad.processing', backend=backend_url, broker=broker_url)


class Parser():
    """
    Instances specify a parser. It allows to find *main files* from  given uploaded 
    and extracted files. Further, allows to run the parser on those 'main files'. 
    """
    def __init__(self, name, main_file_re, main_contents_re):
        self.name = name
        self._main_file_re = re.compile(main_file_re)
        self._main_contents_re = re.compile(main_contents_re)

    def matches(self, upload, filename):
        if self._main_file_re.match(filename):
            try:
                file = upload.open(filename)
                return self._main_contents_re.match(file.read(500))
            finally:
                file.close()


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

@app.task()
def process(upload_id):
    mainfiles = list()
    with files.upload(upload_id) as upload:
        for filename in upload.filelist:
            for parser in parsers:
                if parser.matches(upload, filename):
                    mainfiles.append((filename, parser.name))

    return group([parse.s(mainfile, parser) for mainfile, parser in mainfiles]).delay()


@app.task()
def parse(mainfile, parser):
    return 'parsed %s with %s' % (mainfile, parser)


if __name__ == '__main__':
    print(~process.s('examples_vasp.zip'))
