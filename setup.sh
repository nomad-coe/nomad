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

./generate_gui_artifacts.sh
