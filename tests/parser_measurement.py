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
if __name__ == '__main__':
    import sys
    import logging
    import time
    import json

    from nomad import utils, datamodel
    from nomad.parsing.parsers import parser_dict
    from nomad.cli.parse import normalize_all
    from nomad.metainfo.legacy import LegacyMetainfoEnvironment
    from nomad.parsing.legacy import Backend

    mainfile_path = sys.argv[1]
    utils.set_console_log_level(logging.CRITICAL)

    archive = datamodel.EntryArchive()

    def backend_factory(env, logger):
        return Backend(LegacyMetainfoEnvironment(env), entry_archive=archive, logger=logger)

    logger = utils.get_logger(__name__)
    parser = parser_dict['parsers/vasp']
    setattr(parser, 'backend_factory', backend_factory)

    def run_benchmark():
        for _ in range(0, 10):
            parser.parse(mainfile_path, logger=logger)
            normalize_all(archive)
            with open('/dev/null', 'wt') as f:
                json.dump(archive.m_to_dict(), f)

    start = time.time()
    run_benchmark()
    print(time.time() - start)
