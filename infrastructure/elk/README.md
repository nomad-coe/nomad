This image is based on the populer elk-stack docker image:
[github](https://github.com/spujadas/elk-docker),
[readthedocs](http://elk-docker.readthedocs.io/).

## Changes
- disabled ssl for beats communication to logstash server
- added tcp input
- simplified elastic search output (don't now how to set metric and other vars yet :-()