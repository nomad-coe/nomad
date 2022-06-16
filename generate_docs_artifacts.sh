#!/bin/sh
rm -rf docs/assets/nomad-oasis*
cd ops/docker-compose
zip -r ../../docs/assets/nomad-oasis.zip nomad-oasis -x "**-with-keycloak**" -x "**.empty"
cd ../..