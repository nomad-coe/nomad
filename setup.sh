#!/bin/bash

set -e

pip install --upgrade pip

git submodule sync
sleep 5

# Install nomad
pip install -e .[all]

# Install sub-modules
git submodule update --init --jobs=4
./dependencies.sh -e

./generate_gui_artifacts.sh
