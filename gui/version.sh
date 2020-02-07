#!/bin/sh
commit=`git rev-parse --short --verify HEAD`
sed -i -e "s/nomad-gui-commit-placeholder/$commit/g" package.json
rm -f package.json-e
