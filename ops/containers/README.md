## Containers for additional services

Some third party services cannot be used with ready made public dockerhub images. For those
we provide our own Dockerfiles that extend public base images. These are largely due to
NOMAD specific configuration that we want to apply to the the base images.

### Keycloak (optional)

A custom version of [jboss/keycloak](https://hub.docker.com/r/jboss/keycloak/)
- added bcrypt password hashing: [https://github.com/leroyguillaume/keycloak-bcrypt](https://github.com/leroyguillaume/keycloak-bcrypt)
- create admin user if not there
- import test realm on first use
- use H2 database
- change config to allow reverse proxy under custom prefix


### ELK (optional)

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

### CI runner (optional)

This is the immage that this project uses for its gitlab-ci runner. To build an
push it, you have to log into the project's registry (see [gitlab docs](https://docs.gitlab.com/ee/user/packages/container_registry/)) and do
```
cd ci-runner
docker build -t gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-fair/ci-runner .
docker push gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-fair/ci-runner
```

This image allows to bash, git, docker, docker-compose, k8s, and helm3.