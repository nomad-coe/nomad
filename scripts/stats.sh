#!/bin/sh

set -e

working_dir=$(pwd)
project_dir=$(dirname $(dirname $(realpath $0)))

cd $project_dir

echo "LOC with pygount (pip install pygount)"

echo "backend:       $(pygount nomad/ -s py | awk 'NF {print $1}' | paste -sd+ - | bc)"
echo "backend tests: $(pygount tests/ -s py | awk 'NF {print $1}' | paste -sd+ - | bc)"
echo "frontend:      $(pygount gui/src -s js | awk 'NF {print $1}' | paste -sd+ - | bc)"
