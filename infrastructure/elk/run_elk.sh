#!/bin/bash
docker build -t nomad_elk .
docker run -v nomad_elk:/var/lib/elasticsearch -p 15601:5601 -p 15000:5000 -p 15044:5044 --restart=always nomad_elk