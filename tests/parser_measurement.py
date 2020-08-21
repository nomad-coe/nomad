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
