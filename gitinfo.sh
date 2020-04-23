#!/bin/sh
# python/backend
echo log, ref, version, commit = \"$(git log -1 --oneline)\", \"$(git describe --all)\", \"$(git describe --tags)\", \"$(git rev-parse --verify --short HEAD)\" > nomad/gitinfo.py

# gui
commit=`git rev-parse --short --verify HEAD`
sed -i -e "s/nomad-gui-commit-placeholder/$commit/g" gui/package.json
rm -f gui/package.json-e

python nomad/search.py > gui/src/searchQuantities.json