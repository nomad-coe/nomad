from nomad.dependencies import parser_dict


def test_vasp_parser():
    vasp_parser = parser_dict['parsers/vasp']
    vasp_parser.run('.dependencies/parsers/vasp/test/examples/xml/perovskite.xml')
