import sys
from nomad.client import parse, normalize_all

# match and run the parser
archives = parse(sys.argv[1])
for archive in archives:
    # run all normalizers
    normalize_all(archive)

    # get the 'main section' section_run as a metainfo object
    section_run = archive.run[0]

    # get the same data as JSON serializable Python dict
    python_dict = section_run.m_to_dict()
