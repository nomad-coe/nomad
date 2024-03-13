#!/bin/bash

set -e
set -x # echo on

working_dir=$(pwd)
project_dir=$(dirname $(dirname $(realpath $0)))

cd $project_dir

mkdir tmp

python ops/kubernetes/nomad/updatevalues.py > tmp/helm-values.yaml
diff ops/kubernetes/nomad/values.yaml tmp/helm-values.yaml

# cleanup
rm -rf tmp
