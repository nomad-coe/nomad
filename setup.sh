#!/bin/sh

set -e

git submodule sync
git submodule update --init --jobs=4
./dependencies.sh -e
pip install -e .