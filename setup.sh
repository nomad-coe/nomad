#!/bin/bash

set -e

pip install --upgrade pip

git submodule sync
sleep 5

# Install sub-modules
git submodule update --init --jobs=4
./dependencies.sh -e

# Install nomad
pip install -e .[all]

# Generate GUI artifacts
nomad dev metainfo > gui/src/metainfo.json
nomad dev search-quantities > gui/src/searchQuantities.json
nomad dev units > gui/src/units.js

