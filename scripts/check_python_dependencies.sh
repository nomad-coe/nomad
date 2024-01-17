#!/bin/sh

set -e

working_dir=$(pwd)
project_dir=$(dirname $(dirname $(realpath $0)))

cd $project_dir

set -x # echo on

pip-compile --resolver=backtracking --quiet --annotation-style=line \
    --extra=infrastructure --extra=parsing \
    --output-file=requirements.txt.tmp \
    dependencies/nomad-dos-fingerprints/pyproject.toml \
    dependencies/parsers/atomistic/pyproject.toml \
    dependencies/parsers/database/pyproject.toml \
    dependencies/parsers/eelsdb/pyproject.toml \
    dependencies/parsers/electronic/pyproject.toml \
    dependencies/parsers/nexus/pyproject.toml \
    dependencies/parsers/workflow/pyproject.toml pyproject.toml

diff requirements.txt requirements.txt.tmp


pip-compile --resolver=backtracking --quiet --annotation-style=line \
    --extra=dev --extra=infrastructure --extra=parsing \
    --output-file=requirements-dev.txt.tmp \
    requirements.txt \
    dependencies/nomad-dos-fingerprints/pyproject.toml \
    dependencies/parsers/atomistic/pyproject.toml \
    dependencies/parsers/database/pyproject.toml \
    dependencies/parsers/eelsdb/pyproject.toml \
    dependencies/parsers/electronic/pyproject.toml \
    dependencies/parsers/nexus/pyproject.toml \
    dependencies/parsers/workflow/pyproject.toml \
    pyproject.toml

diff requirements-dev.txt requirements-dev.txt.tmp

# cleanup
rm requirements.txt.tmp requirements-dev.txt.tmp