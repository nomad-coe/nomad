#!/bin/sh
set -e

coe_config_args="-e REPO_DB_JDBC_URL=jdbc:postgresql://postgres:5432/nomad -e REPO_ELASTIC_URL=elasticsearch://elastic:9200 --network nomad_default"

nomad="nomad -p 8000 -h localhost"
repo_tool="docker run $coe_config_args gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-fair/coe-repotool"

echo "reset nomad"
$nomad reset
curl -XDELETE localhost:9200/repo_index
curl -XDELETE localhost:9200/repo_topics

echo "import example calculations"
$nomad upload --unstage tests/data/proc/examples_vasp.zip

echo "create a new index with coe repoTool"
$repo_tool newIndex --indexName=repo_index --indexNameTopics=repo_topics

# try to search for new calculations
curl http://localhost:8111/repo/search/calculation_groups_oldformat?query=repository_program_name%3DVASP