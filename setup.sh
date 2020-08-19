#!/bin/bash

set -e

pip install --upgrade pip

git submodule sync
git submodule update --init --jobs=4
./dependencies.sh -e
pip install -e .[all]

nomad dev metainfo > gui/src/metainfo.json
nomad dev search-quantities > gui/src/search-quantities.json
