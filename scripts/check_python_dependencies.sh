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

# backup
cp requirements.txt requirements.txt.tmp
cp requirements-dev.txt requirements-dev.txt.tmp

pip-compile --resolver=backtracking --annotation-style=line \
    --extra=infrastructure --extra=parsing \
    --output-file=requirements.txt \
    --pip-args="--prefer-binary" \
    dependencies/nomad-dos-fingerprints/pyproject.toml \
    dependencies/parsers/eelsdb/pyproject.toml \
    pyproject.toml

diff requirements.txt.tmp requirements.txt

pip-compile --resolver=backtracking --annotation-style=line \
    --extra=dev --extra=infrastructure --extra=parsing \
    --output-file=requirements-dev.txt \
    --pip-args="--prefer-binary" \
    requirements.txt \
    pyproject.toml

diff requirements-dev.txt.tmp requirements-dev.txt

# cleanup
mv requirements.txt.tmp requirements.txt
mv requirements-dev.txt.tmp requirements-dev.txt
