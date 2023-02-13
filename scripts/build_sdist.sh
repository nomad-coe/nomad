#!/bin/sh
set -e

working_dir=$(pwd)
project_dir=$(dirname $(dirname $(realpath $0)))

cd $project_dir

# Clean up any previous build
rm -rf nomad/app/static/docs
rm -rf nomad/app/static/gui
rm -rf site nomad_lab.egg-info dist build

mkdocs build
mkdir -p nomad/app/static/docs
cp -r site/* nomad/app/static/docs

cd gui
yarn --network-timeout 1200000
yarn run build
cd ..
mkdir -p nomad/app/static/gui
cp -r gui/build/* nomad/app/static/gui

python -m build --sdist
