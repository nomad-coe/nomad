#!/bin/sh
set -e

nomad="docker run gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-fair:test nomad -p 8000 -h api"
repo_tool="docker run gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-fair/coe-repotool:latest"
repo_webservice="docker run gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-fair/coe-repowebservice:latest"

# import example calculations
$nomad upload tests/data/proc/examles_vasp.zip

# create a new index with coe repoTool
$repo_tool newIndex --indexName=repo_index --indexNameTopics=repo_topics

# start coe repoServer
# docker run nomad/coe-repowebservice:latest

# try to search for new calculations
# curl localhost:8111/...