# Development Setup

## Introduction
The nomad infrastructure consists of a series of nomad and 3rd party services:
- nomad worker (python): task worker that will do the processing
- nomad app (python): the nomad app and it's REST APIs
- nomad gui: a small server serving the web-based react gui
- proxy: an nginx server that reverse proxyies all services under one port
- elastic search: nomad's search and analytics engine
- mongodb: used to store processing state
- rabbitmq: a task queue used to distribute work in a cluster

All 3rd party services should be run via *docker-compose* (see blow). The
nomad python  services can be run with python to develop them.
The gui can be run with a development server via yarn.

Below you will find information on how to install all python dependencies and code
manually. How to use *docker*/*docker-compose*. How run 3rd-party services with *docker-compose*.

Keep in mind the *docker-compose* configures all services in a way that mirror
the configuration of the python code in `nomad/config.py` and the gui config in
`gui/.env.development`.

To learn about how to run everything in docker, e.g. to operate a NOMAD OASIS in
production, go (here)(/app/docs/ops.html).

## Install python code and dependencies

### Cloning and development tools
If not already done, you should clone nomad and create a python virtual environment.

To clone the repository:
```
git clone git@gitlab.mpcdf.mpg.de:nomad-lab/nomad-FAIR.git
cd nomad-FAIR
```

### C libs

Even though the NOMAD infrastructure is written in python, there is a C library
required by one of our pyhton dependencies.

#### libmagic

Libmagic allows to determine the MIME type of files. It should be installed on most
unix/linux systems. It can be installed on MacOS with homebrew:

```
brew install libmagic
```

### Virtual environment

#### pyenv
The nomad code currently targets python 3.6. If you host machine has 3.7 or later installed,
you can use [pyenv](https://github.com/pyenv/pyenv) to use python 3.6 in parallel.
To use 3.7 there is a slight issue about the `enum34` which fails the compilation of the
`mdtraj` and `mdanalysis` packages. A possible work arround is to uninstall and tham re-install
`enum34` once the other packages are installed.

#### virtualenv
We strongly recommend to use *virtualenv* to create a virtual environment. It will allow you
to keep nomad and its dependencies separate from your system's python installation.
Make sure to base the virtual environment on Python 3.
To install *virtualenv*, create an environment and activate the environment use:
```
pip install virtualenv
virtualenv -p `which python3` .pyenv
source .pyenv/bin/activate
```

#### Conda
If you are a conda user, there is an equivalent, but you have to install pip and the
right python version while creating the environment.
```
conda create --name nomad_env pip python=3.6
conda activate nomad_env
```

To install libmagick for conda, you can use (other channels might also work):
```
conda -c conda-forge install --name nomad_env libmagic
```

The next steps can be done using the `setup.sh` script. If you prefere to understand all
the steps and run them manually, read on:

### Get all the submodules
We use git submodules to retrieve all the other NOMAD repositories, mainly parsers.

```
git submodules update --init
```

### Install python dependencies
We use *pip* to manage required python packages.
```
pip install -r requirements.txt
```

### Install NOMAD-coe dependencies.
Nomad is based on python modules from the NOMAD-coe project.
This includes parsers, python-common and the meta-info. These modules are maintained as
their own GITLab/git repositories. To clone and initialize them run:

```
git submodules update --init
```

All requirements for these submodules need to be installed and they need to be installed
themselves as python modules. Run the `dependencies.sh` script that will install
everything into your virtual environment:
```
./dependencies.sh -e
```

The `-e` option will install the NOMAD-coe dependencies with symbolic links allowing you
to change the downloaded dependency code without having to reinstall after.

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
the *worker* (does the actual processing), and the *app*. Those services can be
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

### Run necessary 3-rd party services with docker-compose

You can run all containers with:
```
cd ops/docker-compose/nomad
docker-compose -f docker-compose.yml -f docker-compose.override.yml up -d mongo elastic rabbitmq
```

To shut down everything, just `ctrl-c` the running output. If you started everything
in *deamon* mode (`-d`) use:
```
docker-compose down
```

Usually these services only used by the nomad containers, but sometimes you also
need to check something or do some manual steps.

The *docker-compose* can be overriden with additional seetings. See documentation section on
operating NOMAD for more details. The override `docker-compose.override.yml` will
expose all database ports to the hostmachine and should be used in development. To use
it run docker-compose with `-f docker-compose.yml -f docker-compose.override.yml`.

### ELK (elastic stack)

If you run the ELK stack (and enable logstash in nomad/config.py),
you can reach the Kibana with [localhost:5601](http://localhost:5601).
The index prefix for logs is `logstash-`. The ELK is only available with the
`docker-compose.dev-elk.yml` override.

### mongodb and elastic search

You can access mongodb and elastic search via your preferred tools. Just make sure
to use the right ports (see above).

## Run nomad services

### API and worker

To simply run a worker with the installed nomad cli, do (from the root)
```
nomad admin run worker
```

To run it directly with celery, do (from the root)
```
celery -A nomad.processing worker -l info
```

You can also run worker and app together:
```
nomad admin run appworker
```

### GUI
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
docker-compose up -d elastic rabbitmq
cd ../..
pytest -svx tests
```

We use pylint, pycodestyle, and mypy to ensure code quality. To run those:
```
nomad dev qa --skip-test
```

To run all tests and code qa:
```
nomad dev qa
```

This mimiques the tests and checks that the GitLab CI/CD will perform.


## Setup your (I)DE

The documentation section on development guidelines details how the code is organized,
tested, formatted, and documented. To help you meet these guidelines, we recomment to
use a proper IDE for development and ditch any VIM/Emacs (mal-)practices.

### Visual Studio Code

Here are some VSCode settings that will enable features for linting, some auto formating,
line size ruler, etc.
```json
{
    "python.venvPath": "${workspaceFolder}/.pyenv",
    "python.pythonPath": "${workspaceFolder}/.pyenv/bin/python",
    "git.ignoreLimitWarning": true,
    "editor.rulers": [90],
    "editor.renderWhitespace": "all",
    "editor.tabSize": 4,
    "[javascript]": {
        "editor.tabSize": 2
    },
    "files.trimTrailingWhitespace": true,
    "git.enableSmartCommit": true,
    "eslint.autoFixOnSave": true,
    "python.linting.pylintArgs": [
        "--load-plugins=pylint_mongoengine",
    ],
    "python.linting.pep8Path": "pycodestyle",
    "python.linting.pep8Enabled": true,
    "python.linting.pep8Args": ["--ignore=E501,E701"],
    "python.linting.mypyEnabled": true,
    "python.linting.mypyArgs": [
        "--ignore-missing-imports",
        "--follow-imports=silent",
        "--no-strict-optional"
    ],
    "workbench.colorCustomizations": {
        "editorError.foreground": "#FF2222",
        "editorOverviewRuler.errorForeground": "#FF2222",
        "editorWarning.foreground": "#FF5500",
        "editorOverviewRuler.warningForeground": "#FF5500",
        "activityBar.background": "#4D2111",
        "titleBar.activeBackground": "#6B2E18",
        "titleBar.activeForeground": "#FDF9F7"
    },
    "files.watcherExclude": {
        "**/.git/objects/**": true,
        "**/.git/subtree-cache/**": true,
        "**/node_modules/*/**": true,
        "**/.pyenv/*/**": true,
        "**/__pycache__/*/**": true,
        "**/.mypy_cache/*/**": true,
        "**/.volumes/*/**": true,
        "**/docs/.build/*/**": true
    }
}
```

Here are some example launch configs for VSCode:

```json
{
  "version": "0.2.0",
  "configurations": [
    {
      "type": "chrome",
      "request": "launch",
      "name": "Launch Chrome against localhost",
      "url": "http://localhost:3000",
      "webRoot": "${workspaceFolder}/gui"
    },
    {
      "name": "Python: API Flask (0.11.x or later)",
      "type": "python",
      "request": "launch",
      "module": "flask",
      "env": {
        "FLASK_APP": "nomad/app/__init__.py"
      },
      "args": [
        "run",
        "--port",
        "8000",
        "--no-debugger",
        "--no-reload"
      ]
    },
    {
      "name": "Python: some test",
      "type": "python",
      "request": "launch",
      "cwd": "${workspaceFolder}",
      "program": "${workspaceFolder}/.pyenv/bin/pytest",
      "args": [
        "-sv",
        "tests/test_cli.py::TestClient::test_mirror"
      ]
    },
    {
      "name": "Python: Current File",
      "type": "python",
      "request": "launch",
      "program": "${file}"
    },
    {
      "name": "Python: Attach",
      "type": "python",
      "request": "attach",
      "localRoot": "${workspaceFolder}",
      "remoteRoot": "${workspaceFolder}",
      "port": 3000,
      "secret": "my_secret",
      "host": "localhost"
    }
  ]
}
```