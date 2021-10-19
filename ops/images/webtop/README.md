This is a test image for north tools authorization.

To build:
```
docker build -t gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-remote-tools-hub/webtop .
```

To push:
```
docker login gitlab-registry.mpcdf.mpg.de
docker push gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-remote-tools-hub/webtop
```

To use it, you have to adopt `/nomad/jupyterhub_config.py` to use the
`gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-remote-tools-hub/webtop` as an image.