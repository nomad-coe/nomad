#!/bin/sh
set -e

working_dir=$(pwd)
project_dir=$(dirname $(dirname $(realpath $0)))

cd $project_dir

# Clean up any previous build
rm -rf site nomad/app/static/docs
mkdocs build
mkdir -p nomad/app/static/docs
cp -r site/* nomad/app/static/docs/

rm -rf nomad/app/static/docs/gui
cd gui
yarn --network-timeout 1200000
yarn run build
rm build/env.js
mkdir -p ../nomad/app/static/gui
cp -r build/* ../nomad/app/static/gui
cd ..

rm -rf site nomad_lab.egg-info dist build
python -m build