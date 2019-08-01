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
from sys import stdout

from nomad import processing as proc, search, datamodel, infrastructure

from nomad.cli.cli import cli


@cli.group(help='''The nomad admin commands to do nasty stuff directly on the databases.
                     Remember: With great power comes great responsibility!''')
@click.pass_context
def admin(ctx):
    pass


@admin.command(help='(Re-)index all calcs.')
def index():
    infrastructure.setup_logging()
    infrastructure.setup_mongo()
    infrastructure.setup_elastic()
    infrastructure.setup_repository_db()

    all_calcs = proc.Calc.objects().count()
    print('indexing %d ...' % all_calcs)

    def calc_generator():
        count = 0

        for calc in proc.Calc.objects():
            count += 1
            if count % 1000 == 0:
                stdout.write('\r   indexed %10d of %10d calcs' % (count, all_calcs))
                stdout.flush()
            yield datamodel.CalcWithMetadata(**calc.metadata)

    failed = search.index_all(calc_generator())
    print('')
    print('indexing completed, %d failed entries' % failed)
