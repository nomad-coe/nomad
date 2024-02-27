#!/bin/bash

interval=8

while true; do
  echo "Checking if the cluster is green in $interval seconds..."
  sleep $interval
  result=$(curl -s http://elastic:9200/_cat/health | awk '{print $4}')
  if [ "$result" == "green" ]; then
    exit 0
  fi
  # double the interval if the cluster is not green
  interval=$((interval * 2))
  if [ $interval -gt 120 ]; then
    echo "Cluster is not green after 2 minutes."
    exit 1
  fi
done
