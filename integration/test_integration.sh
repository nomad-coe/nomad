#!/bin/sh
set -e

coe_config_args="-e \"REPO_DB_JDBC_URL=jdbc:postgresql://postgres:5432/nomad\" -e \"REPO_ELASTIC_URL=elasticsearch://elastic:9200\""

nomad="docker run gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-fair:test python -m nomad.client -p 8000 -h api"
repo_tool="docker $coe_config_args run gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-fair/coe-repotool:latest"
repo_webservice="docker $coe_config_args run gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-fair/coe-repowebservice:latest"

# import example calculations
$nomad upload --unstage /app/tests/data/proc/examples_vasp.zip

# create a new index with coe repoTool
$repo_tool newIndex --indexName=repo_index --indexNameTopics=repo_topics

# start coe repoServer
# docker run nomad/coe-repowebservice:latest

# try to search for new calculations
# curl localhost:8111/...