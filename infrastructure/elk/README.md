# ELK

This image is based on the populer elk-stack docker image:
[github](https://github.com/spujadas/elk-docker),
[readthedocs](http://elk-docker.readthedocs.io/).

## Changes
- disabled ssl for beats communication to logstash server
- added tcp input
- simplified elastic search output (don't now how to set metric and other vars yet :-()
- added kibana.yml::server.basePath="/nomad/kibana"


## Usage
You can run this image outside the usual docker-compose.

To use this image with reverse proxy in nginx, use:

```
location ~ ^/nomad/kibana/(.*)$ {
    proxy_pass http://130.183.207.116:15601/$1;
    proxy_set_header Host $host;
    proxy_set_header X-Real-IP $remote_addr;
}
```