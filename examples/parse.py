import sys

from nomad import utils
from nomad.cli.parse import parse


utils.configure_logging()
parse(sys.argv[1])
