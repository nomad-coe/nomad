#!/bin/sh
rm -rf docs/assets/nomad-oasis*
cd ops/docker-compose
zip -r ../../docs/assets/nomad-oasis.zip nomad-oasis -x "**/.gitignore"
zip -r ../../docs/assets/nomad-oasis-with-keycloak.zip nomad-oasis-with-keycloak -x "**/.gitignore"
cd ../..