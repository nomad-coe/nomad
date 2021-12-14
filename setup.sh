#!/bin/bash

set -e

pip install --upgrade pip

git submodule sync --recursive
sleep 5

# Install sub-modules
git submodule update --init --recursive --jobs=4
./dependencies.sh -e

# Install nomad
pip install -e .[all]

./generate_gui_artifacts.sh
