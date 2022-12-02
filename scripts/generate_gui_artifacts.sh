#!/bin/sh

set -e

working_dir=$(pwd)
project_dir=$(dirname $(dirname $(realpath $0)))

cd $project_dir

python -m nomad.cli dev gui-artifacts --output-directory gui/src

NOMAD_CONFIG=gui/tests/nomad.yaml python -m nomad.cli dev gui-config >gui/tests/env.js
