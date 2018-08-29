#!/bin/bash
docker build -t nomadxt_elk .
docker run -v nomadxt_elk:/var/lib/elasticsearch -p 15601:5601 -p 15000:5000 -p 15044:5044 --restart=always nomadxt_elk