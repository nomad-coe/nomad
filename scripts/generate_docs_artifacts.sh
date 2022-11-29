#!/bin/sh

set -e

working_dir=$(pwd)
project_dir=$(dirname $(dirname $(realpath $0)))

cd $project_dir/ops/docker-compose

rm -rf $project_dir/docs/assets/nomad-oasis*
zip -r $project_dir/docs/assets/nomad-oasis.zip nomad-oasis -x "**/.gitignore"
zip -r $project_dir/docs/assets/nomad-oasis-with-keycloak.zip nomad-oasis-with-keycloak -x "**/.gitignore"
