#!/bin/bash

set -e

pip install --upgrade pip

git submodule sync --recursive
sleep 5

# Install nomad
pip install -e .[all]

# Install sub-modules
git submodule update --init --recursive --jobs=4
./dependencies.sh -e

./generate_gui_artifacts.sh

rm -rf docs/build && mkdocs build && mv site docs/build
