#!/bin/sh

set -e

git submodule sync
git submodule update --init
pip install -r requirements.txt
./dependencies.sh -e
pip install -e .