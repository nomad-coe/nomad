# Development Setup

## Introduction
The nomad infrastructure consists of a series of nomad and 3rd party services:
- nomad worker (python): task worker that will do the processing
- nomad api (python): the nomad REST API
- nomad gui: a small server serving the web-based react gui
- proxy: an nginx server that reverse proxyies all services under one port
- elastic search: nomad's search and analytics engine
- mongodb: used to store processing state
- rabbitmq: a task queue used to distribute work in a cluster

All 3rd party services should be run via *docker-compose* (see blow). The
nomad python  services can also be run via *docker-compose* or manually started with python.
The gui can be run manually with a development server via yarn, or with
*docker-compose*

Below you will find information on how to install all python dependencies and code
manually. How to use *docker*/*docker-compose*. How run services with *docker-compose*
or manually.

Keep in mind the *docker-compose* configures all services in a way that mirror
the configuration of the python code in `nomad/config.py` and the gui config in
`gui/.env.development`.

## Install python code and dependencies

### Cloning and development tools
If not already done, you should clone nomad and create a python virtual environment.

To clone the repository:
```
git clone git@gitlab.mpcdf.mpg.de:nomad-lab/nomad-FAIR.git
cd nomad-FAIR
```

The nomad code currently targets python 3.6. If you host machine has 3.7 or later installed,
you can use [pyenv](https://github.com/pyenv/pyenv) to use python 3.6 in parallel.
While in principle everything should be compatable with 3.7 and later there have been
issues with some dependencies and requirements not being compatible with 3.7

We strongly recommend to use *virtualenv* to create a virtual environment. It will allow you
to keep nomad and its dependencies separate from your system's python installation.
Make sure to base the virtual environment on Python 3.
To install *virtualenv*, create an environment and activate the environment use:
```
pip install virtualenv
virtualenv -p `which python3` .pyenv
source .pyenv/bin/activate
```

We use *pip* to manage required python packages.
```
pip install -r requirements.txt
```

### Install NOMAD-coe dependencies.
Nomad is based on python modules from the NOMAD-coe project.
This includes parsers, normalizers, python-common and the meta-info.
Those dependencies are managed and configured via python in
`nomad/dependencies.py`. This gives us more flexibility in interacting with
different parser, normalizer versions from within the running nomad infrastructure.

To run the dependencies script and install all dependencies into your environment:
```
./dependencies.sh
```
This will checkout the proper version of the respective NOMAD-coe modules, install
further requirements, and install the modules themselves. The `-e` option will install
the NOMAD-coe dependencies with symbolic links allowing you to change the downloaded
dependency code without having to reinstall after.

### Install nomad
Finally, you can add nomad to the environment itself.
```
pip install -e .
```

## Build and run the infrastructure with docker

### Docker and nomad
Nomad depends on a set of databases, search engines, and other services. Those
must run to make use of nomad. We use *docker* and *docker-compose* to create a
unified environment that is easy to build and to run.

You can use *docker* to run all necessary 3rd-party components and run all nomad
services manually from your python environment. Or you can run everything within
docker containers. The former is often preferred during development, since it allows
you change things, debug, and re-run things quickly. The later one brings you
closer to the environment that will be used to run nomad in production.

### Docker images for nomad
There are currently two different images and respectively two different docker files:
`Dockerfile`, and `gui/Dockerfile`.

Nomad comprises currently two services,
the *worker* (does the actual processing), and the *api*. Those services can be
run from one image that have the nomad python code and all dependencies installed. This
is covered by the `Dockerfile`.

The gui is served via containers based on the `gui/Dockerfile` which contains the
react-js frontend code. Before this image can be build, make sure to execute

```
cd gui
./gitinfo.sh
cd ..
```

This allows to gui to present some information about the current git revision without
having to copy the git itself to the docker build context.

The images are build via *docker-compose* and don't have to be created manually.

### Build with docker-compose

We have multiple *docker-compose* files that must be used together.
- `docker-compose.yml` contains the base definitions for all services
- `docker-compose.override.yml` configures services for development (notably builds images for nomad services)
- `docker-compose.prod.yml` configures services for production (notable uses a pre-build image for nomad services that was build during CI/CD)

It is sufficient to use the implicit `docker-compose.yml` only (like in the command below).
The `override` will be used automatically.

There is also an `.env` file. For development you can use `.env_development`:
```
cd ./ops/docker-compose/nomad
ln -s .env_development .env
```

The production `.env` file is stored on our serves and not part of the source code.

Now we can build the *docker-compose* that contains all external services (rabbitmq,
mongo, elastic, elk) and nomad services (worker, api, gui).
```
docker-compose build
```

Docker-compose tries to cache individual building steps. Sometimes this causes
troubles and not everything necessary is build when you changed something. In
this cases use:
```
docker-compose build --no-cache
```

### Run everything with docker-compose

You can run all containers with:
```
docker-compose up
```

To shut down everything, just `ctrl-c` the running output. If you started everything
in *deamon* mode (`-d`) use:
```
docker-compose down
```

### Run containers selectively
The following services/containers are managed via our docker-compose:
- rabbitmq, mongo, elastic, (elk, only for production)
- worker, api
- gui
- proxy

The *proxy* container runs *nginx* based reverse proxies that put all services under
a single port and different paths.

You can also run services selectively, e.g.
```
docker-compose up -d rabbitmq, mongo, elastic, elk
docker-compose up worker
docker-compose up api gui proxy
```

### Configure the containers
The *docker-compose* takes some configuration from the environment. Environment
variables are set in `.env`, which is a link to either `.env_local` (intended to run
nomad for development on a local computer) and `.env_processing` (indented to run
nomad on nomad in *'production'*, currently on the enc pre-processing machine).

You can configure host ports, volume locations for host bindings, and these sort of things.

## Accessing 3'rd party services

Usually these services only used by the nomad containers, but sometimes you also
need to check something or do some manual steps.

The file `ops/docker-compose/nomad/.env` contains variables that control the ports
used to bind internal docker ports to your host machine. These are the ports you
have to use to connect to the respective services.

### ELK (elastic stack)

If you run the ELK stack (and enable logstash in nomad/config.py),
you can reach the Kibana with [localhost:5601](http://localhost:5601).
The index prefix for logs is `logstash-`.

### mongodb and elastic search

You can access mongodb and elastic search via your preferred tools. Just make sure
to use the right ports (see above).


## Run nomad services manually

You can run the worker, api, and gui as part of the docker infrastructure, like
seen above. But, of course there are always reasons to run them manually during
development, like running them in a debugger, profiler, etc.

### Run the nomad worker manually

To simply run a worker with the installed nomad cli, do (from the root)
```
nomad run worker
```

To run it manually with celery, do (from the root)
```
celery -A nomad.processing worker -l info
```
You can use different debug level (e.g. switch `info` to `debug`)

Use watchdog during development to reload the worker on code changes.
Watchdog is part of the requirements-dev.txt. For MacOS (there is currently a bug in watchdog)
uninstall and install this [fixed](https://github.com/gorakhargosh/watchdog/issues/330) version
```
pip uninstall watchdog
pip install git+https://github.com/gorakhargosh/watchdog.git
```

Now use this to auto relead worker:
```
watchmedo auto-restart -d ./nomad -p '*.py' -- celery worker -l info -A nomad.processing
```

### Run the api
Either with docker, or:
```
nomad run api
```

Or manually:
```
python nomad/api.py
```

### Run the gui
When you run the gui on its own (e.g. with react dev server below), you have to have
the API running manually also. This *inside docker* API is configured for ngingx paths
and proxies, which are run by the gui container. But you can run the *production* gui
in docker and the dev server gui in parallel with an API in docker.
Either with docker, or:
```
cd gui
yarn
yarn start
```

## Run the tests
You need to have the infrastructure partially running: elastic, rabbitmq.
The rest should be mocked or provided by the tests. Make sure that you do no run any
worker, as they will fight for tasks in the queue.
```
cd ops/docker-compose
docker-compose up -d elastic rabbitmq postgres
cd ../..
pytest -svx tests
```

We use pylint, pycodestyle, and mypy to ensure code quality. To run those:
```
nomad qa --skip-test
```

To run all tests and code qa:
```
nomad qa
```

This mimiques the tests and checks that the GitLab CI/CD will perform.
