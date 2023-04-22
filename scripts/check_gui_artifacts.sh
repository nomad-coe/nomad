#!/bin/bash

set -e

working_dir=$(pwd)
project_dir=$(dirname $(dirname $(realpath $0)))

cd $project_dir

set -x # echo on

mkdir tmp

NOMAD_CONFIG=gui/tests/nomad.yaml python -m nomad.cli dev gui-artifacts > tmp/artifacts.js
NOMAD_CONFIG=gui/tests/nomad.yaml python -m nomad.cli dev gui-config > tmp/env.js

diff gui/tests/artifacts.js tmp/artifacts.js
diff gui/tests/env.js tmp/env.js

# cleanup
rm -rf tmp


