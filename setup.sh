#!/bin/bash

set -e

pip install --upgrade pip

git submodule sync
sleep 5
git submodule update --init --jobs=4
./dependencies.sh -e
pip install -e .[all]

./generate_gui_artifacts.sh