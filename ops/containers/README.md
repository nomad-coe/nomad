## Containers

### ELK

This image is based on the popular elk-stack docker image:
[github](https://github.com/spujadas/elk-docker),
[readthedocs](http://elk-docker.readthedocs.io/).

Changes
- disabled ssl for beats communication to logstash server
- added tcp input
- simplified elastic search output (don't now how to set metric and other vars yet :-()
- added kibana.yml::server.basePath="/nomad/kibana"

The file `elk/kibana_objects.json` contains an export of nomad specific searches,
visualizations, and dashboard.