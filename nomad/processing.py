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
import nomad.config as config
import nomad.files as files

broker_url = 'pyamqp://%s:%s@localhost//' % (config.rabbitmq.user, config.rabbitmq.password)
backend_url = 'rpc://localhost'
app = Celery('nomad.processing', backend=backend_url, broker=broker_url)


@app.task()
def process(upload):
  mainfiles = [('a', 'pa'), ('b', 'pb'), ('c', 'pc')]
  parsers = group([parse.s(mainfile, parser) for mainfile, parser in mainfiles]).delay()
  return parsers

@app.task()
def parse(mainfile, parser):
  return 'parsed %s with %s' % (mainfile, parser)

if __name__ == '__main__':
  print(~process.s('test'))
