#!/bin/bash

set -e

working_dir=$(pwd)
project_dir=$(dirname $(dirname $(realpath $0)))

cd $project_dir

# Initialise all of the submodules
git submodule update --init --recursive

# Clean up any previous build
rm -rf nomad/app/static/docs
rm -rf nomad/app/static/gui
rm -rf site

# Install nomad
pip install -e ".[infrastructure,parsing,dev]" -c requirements-dev.txt

# Install "default" plugins. TODO: This can be removed once we have proper
# distributions projects.
pip install -r default_plugins.txt -c requirements-dev.txt

# Build documentation
sh scripts/generate_docs_artifacts.sh
mkdocs build
mkdir -p nomad/app/static/docs
cp -r site/* nomad/app/static/docs

# Generate .env file for the GUI
python -m nomad.cli dev gui-env > gui/.env.development
