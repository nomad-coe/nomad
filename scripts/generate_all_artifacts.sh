#!/bin/sh

set -e

working_dir=$(pwd)
project_dir=$(dirname $(dirname $(realpath $0)))

cd $project_dir/scripts

./generate_docs_artifacts.sh
./generate_gui_test_artifacts.sh
./generate_example_uploads.sh

