#!/bin/bash

set -e

working_dir=$(pwd)
project_dir=$(dirname $(dirname $(realpath $0)))

cd $project_dir

# Install nomad
pip install --prefer-binary -r requirements-dev.txt

# Build documentation
mkdocs build
mkdir -p nomad/app/static/docs
cp -r site/* nomad/app/static/docs/

# Build and copy gui code
cd gui
yarn --network-timeout 1200000
yarn run build
cd ..

cp gui/build nomad/app/static/gui
rm nomad/app/static/gui/env.js
