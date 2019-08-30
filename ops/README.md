## Overview

Read the [introduction](./introduction.html) and [setup](./setup.html) for input on
the different nomad services. This is about how to deploy and operate these services.

The databases and other external services can be run from
normal public dockerhub images or you have to setup your own private registry with
respective images.

The NOMAD specific images are provide by our
[gitlab container registry](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/container_registry).
To access these container you have to login without MPCDF gitlab account:

```
docker login gitlab-registry.mpcdf.mpg.de/nomad-lab
```

There are two basic options to run all these containers: docker-compose and kubernetes.