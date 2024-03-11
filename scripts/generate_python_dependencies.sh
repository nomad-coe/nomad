#!/bin/sh

set -e

working_dir=$(pwd)
project_dir=$(dirname $(dirname $(realpath $0)))

cd $project_dir

pip-compile --resolver=backtracking --annotation-style=line \
    --extra=infrastructure --extra=parsing \
    --output-file=requirements.txt \
    dependencies/nomad-dos-fingerprints/pyproject.toml \
    dependencies/parsers/atomistic/pyproject.toml \
    dependencies/parsers/database/pyproject.toml \
    dependencies/parsers/eelsdb/pyproject.toml \
    dependencies/parsers/electronic/pyproject.toml \
    dependencies/parsers/nexus/pyproject.toml \
    dependencies/parsers/workflow/pyproject.toml pyproject.toml


pip-compile --resolver=backtracking --annotation-style=line \
    --extra=dev --extra=infrastructure --extra=parsing \
    --output-file=requirements-dev.txt \
    requirements.txt \
    pyproject.toml

