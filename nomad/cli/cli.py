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
import logging
import os

from nomad import config as nomad_config, infrastructure


@click.group(help='''This is the entry point to nomad\'s command line interface CLI.
                     It uses a sub-command structure similar to the git command.''')
@click.option('-v', '--verbose', help='sets log level to info', is_flag=True)
@click.option('--debug', help='sets log level to debug', is_flag=True)
@click.option('--config', help='the config file to use')
@click.pass_context
def cli(ctx, verbose: bool, debug: bool, config: str):
    if config is not None:
        nomad_config.load_config(config_file=config)

    if debug:
        nomad_config.console_log_level = logging.DEBUG
    elif verbose:
        nomad_config.console_log_level = logging.INFO
    else:
        nomad_config.console_log_level = logging.WARNING

    nomad_config.service = os.environ.get('NOMAD_SERVICE', 'admin')
    infrastructure.setup_logging()
