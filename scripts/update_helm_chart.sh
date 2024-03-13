#!/bin/sh

set -e

working_dir=$(pwd)
project_dir=$(dirname $(dirname $(realpath $0)))

cd $project_dir

python ops/kubernetes/nomad/updatevalues.py ops/kubernetes/nomad/values.yaml
