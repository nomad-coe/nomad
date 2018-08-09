import unittest
from unittest import TestCase
import logging

from nomad.dependencies import parser_dict


class DependenciesTests(TestCase):

    def test_vasp_parser(self):
        vasp_parser = parser_dict['parsers/vasp']
        vasp_parser.run('.dependencies/parsers/vasp/test/examples/xml/perovskite.xml')

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    unittest.main()
