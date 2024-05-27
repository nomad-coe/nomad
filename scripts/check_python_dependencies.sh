#!/bin/sh

# For now we exclude the following dependencies, because they import an older version
# of nomad-lab causing inevitable conflicts.
#    dependencies/parsers/atomistic/pyproject.toml \
#    dependencies/parsers/database/pyproject.toml \
#    dependencies/parsers/electronic/pyproject.toml \
#    dependencies/parsers/workflow/pyproject.toml \

set -e

working_dir=$(pwd)
project_dir=$(dirname $(dirname $(realpath $0)))

cd $project_dir

# do no support git install for pynxtools
grep pynxtools pyproject.toml | grep "@git+" 1>&2

# backup
cp requirements.txt requirements.txt.tmp
cp requirements-dev.txt requirements-dev.txt.tmp

uv pip compile -q -U --annotation-style=line \
    --extra=infrastructure --extra=parsing \
    --output-file=requirements.txt \
    dependencies/nomad-dos-fingerprints/pyproject.toml \
    dependencies/parsers/eelsdb/pyproject.toml \
    pyproject.toml

diff requirements.txt.tmp requirements.txt

uv pip compile -q -U --annotation-style=line \
    --extra=dev --extra=infrastructure --extra=parsing \
    --output-file=requirements-dev.txt \
    requirements.txt \
    pyproject.toml

diff requirements-dev.txt.tmp requirements-dev.txt

# cleanup
mv requirements.txt.tmp requirements.txt
mv requirements-dev.txt.tmp requirements-dev.txt
