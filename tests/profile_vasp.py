from nomad.parsing import parser_dict
from nomad import utils


def run_parser(parser_name, mainfile):
    parser = parser_dict[parser_name]
    return parser.run(mainfile, logger=utils.get_logger(__name__))


if __name__ == '__main__':
    run_parser('parsers/vasp', '.dependencies/parsers/vasp/test/examples/xml/perovskite.xml')
    run_parser('parsers/vasp', '.dependencies/parsers/vasp/test/examples/xml/perovskite.xml')
    run_parser('parsers/vasp', '.dependencies/parsers/vasp/test/examples/xml/perovskite.xml')
    run_parser('parsers/vasp', '.dependencies/parsers/vasp/test/examples/xml/perovskite.xml')
    run_parser('parsers/vasp', '.dependencies/parsers/vasp/test/examples/xml/perovskite.xml')
