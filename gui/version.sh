#!/bin/sh
version=`git rev-parse --verify HEAD`
sed -i -e "s/nomad-gui-version-placeholder/$version/g" package.json
rm -f package.json-e
