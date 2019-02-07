#!/bin/sh

set -e

echo "log into docker registry..."
docker login gitlab-registry.mpcdf.mpg.de -u $1 -p $2

echo "building images..."
cd dependencies/nomad-lab-base
sbt "project repoTool" docker
sbt "project repoWebservice" docker

echo "pushing images..."
docker push gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-fair/coe-repotool
docker push gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-fair/coe-repowebservice