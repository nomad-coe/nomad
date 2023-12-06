#!/bin/sh

set -e

working_dir=$(pwd)
project_dir=$(dirname $(dirname $(realpath $0)))

cd $project_dir/ops/docker-compose

rm -rf $project_dir/docs/assets/nomad-oasis*
zip -r $project_dir/docs/assets/nomad-oasis.zip nomad-oasis -x "**/.gitignore"
zip -r $project_dir/docs/assets/nomad-oasis-with-keycloak.zip nomad-oasis-with-keycloak -x "**/.gitignore"

rm -rf nomad-oasis-with-plugins/.volumes nomad-oasis-with-plugins/configs nomad-oasis-with-plugins/docker-compose.yaml nomad-oasis-with-plugins/nomad-*
cp -r nomad-oasis/configs nomad-oasis/docker-compose.yaml nomad-oasis/.volumes nomad-oasis-with-plugins/
git clone https://github.com/nomad-coe/nomad-schema-plugin-example.git ./nomad-oasis-with-plugins/nomad-schema-plugin-example
git clone https://github.com/nomad-coe/nomad-parser-plugin-example.git ./nomad-oasis-with-plugins/nomad-parser-plugin-example
git clone https://github.com/nomad-coe/nomad-normalizer-plugin-example.git ./nomad-oasis-with-plugins/nomad-normalizer-plugin-example
cat nomad-oasis-with-plugins/nomad.yaml >> nomad-oasis-with-plugins/configs/nomad.yaml
zip -r $project_dir/docs/assets/nomad-oasis-with-plugins.zip nomad-oasis-with-plugins -x "**/.gitignore" -x "nomad-oasis-with-plugins/nomad.yaml"