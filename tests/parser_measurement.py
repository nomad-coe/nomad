if __name__ == '__main__':
    import sys
    import logging
    import time
    from nomad import config, utils
    from nomad.parsing import parser_dict
    from nomad.cli.parse import normalize_all
    from nomad.metainfo.legacy import LegacyMetainfoEnvironment
    from nomad.parsing.legacy import Backend

    mainfile_path = sys.argv[1]
    utils.set_console_log_level(logging.CRITICAL)

    def backend_factory(env, logger):
        return Backend(LegacyMetainfoEnvironment(env), logger=logger)

    logger = utils.get_logger(__name__)
    parser = parser_dict['parsers/vasp']
    setattr(parser, 'backend_factory', backend_factory)

    def run_benchmark():
        for _ in range(0, 10):
            backend = parser.run(mainfile_path, logger=logger)

            if not backend.status[0] == 'ParseSuccess':
                logger.error('parsing was not successful', status=backend.status)

            backend.openNonOverlappingSection('section_entry_info')
            backend.addValue('upload_id', config.meta.services.unavailable_value)
            backend.addValue('calc_id', config.meta.services.unavailable_value)
            backend.addValue('calc_hash', "no hash")
            backend.addValue('mainfile', mainfile_path)
            backend.addValue('parser_name', 'parsers/vasp')
            backend.closeNonOverlappingSection('section_entry_info')

            normalize_all(backend)
            with open('/dev/null', 'wt') as f:
                backend.write_json(f, pretty=True)

    start = time.time()
    run_benchmark()
    print(time.time() - start)
