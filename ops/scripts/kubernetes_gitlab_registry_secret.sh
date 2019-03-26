#!/bin/bash
kubectl create secret docker-registry gitlab-mpcdf \
  --docker-server https://gitlab-registry.mpcdf.mpg.de \
  --docker-username $1 \
  --docker-password $2 \
  --docker-email $3
