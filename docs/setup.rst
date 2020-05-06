Development Setup
=================
Introduction
------------
The nomad infrastructure consists of a series of nomad and 3rd party services:

 - worker (python): task worker that will do the processing
 - app (python): the nomad app, it's REST APIs and the GUI
 - elastic: nomad's search and analytics engine
 - mongodb: used to store processing state
 - rabbitmq: a task queue used to distribute work in a cluster
 - keycloak: centralized user management

All 3rd party services should be run via *docker-compose* (see below). The
nomad python services can be run with python to develop them. The gui can be
run with a development server via yarn.

Below you will find information on how to install all python dependencies and
code manually, how to run 3rd-party services with *docker-compose*, and how run
and interact with the services in an development environment. To learn about
how to run everything in a production envrionment under NOMAD OASIS, go to :doc:`oasis`.

Keep in mind the *docker-compose* configures all services in a way that mirror
the configuration of the python code in `nomad/config.py` and the gui config in
`gui/.env.development`.


Install python code and dependencies
------------------------------------
Cloning and development tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If not already done, you should clone nomad and create a python virtual environment.

To clone the repository:

.. code-block:: sh

    git clone git@gitlab.mpcdf.mpg.de:nomad-lab/nomad-FAIR.git
    cd nomad-FAIR

C libs
^^^^^^
Even though the NOMAD infrastructure is written in python, there is a C library
required by one of our pyhton dependencies.

libmagic
~~~~~~~~
Libmagic allows to determine the MIME type of files. It should be installed on most
unix/linux systems. It can be installed on MacOS with homebrew:

.. code-block:: sh

    brew install libmagic

Virtual environment
^^^^^^^^^^^^^^^^^^^
pyenv
~~~~~
The nomad code currently targets python 3.6. If you host machine has 3.7 or later installed,
you can use `pyenv <https://github.com/pyenv/pyenv>`_ to use python 3.6 in parallel.
To use 3.7 there is a slight issue about the :code:`enum34` which fails the compilation of the
:code:`mdtraj` and :code:`mdanalysis` packages. A possible work arround is to uninstall and tham re-install
:code:`enum34` once the other packages are installed.

virtualenv
~~~~~~~~~~
We strongly recommend to use *virtualenv* to create a virtual environment. It will allow you
to keep nomad and its dependencies separate from your system's python installation.
Make sure to base the virtual environment on Python 3.
To install *virtualenv*, create an environment and activate the environment use:

.. code-block:: sh

    pip install virtualenv
    virtualenv -p `which python3` .pyenv
    source .pyenv/bin/activate

Conda
~~~~~
If you are a conda user, there is an equivalent, but you have to install pip and the
right python version while creating the environment.

.. code-block:: sh

    conda create --name nomad_env pip python=3.6
    conda activate nomad_env

To install libmagick for conda, you can use (other channels might also work):

.. code-block:: sh

    conda -c conda-forge install --name nomad_env libmagic

The next steps can be done using the :code:`setup.sh` script. If you prefere to understand all
the steps and run them manually, read on:

Install python dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^
We use *pip* to manage required python packages.

.. code-block:: sh

    pip install -r requirements.txt

Install NOMAD-coe dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Nomad is based on python modules from the NOMAD-coe project.
This includes parsers, python-common and the meta-info. These modules are maintained as
their own GITLab/git repositories. To clone and initialize them run:

.. code-block:: sh

    git submodule update --init

All requirements for these submodules need to be installed and they need to be installed
themselves as python modules. Run the :code:`dependencies.sh` script that will install
everything into your virtual environment:

.. code-block:: sh

    ./dependencies.sh -e

The :code:`-e` option will install the NOMAD-coe dependencies with symbolic links allowing you
to change the downloaded dependency code without having to reinstall after.

Install nomad
^^^^^^^^^^^^^
Finally, you can add nomad to the environment itself.

.. code-block:: sh

    pip install -e .

Build and run the full development infrastructure with docker-compose
---------------------------------------------------------------------
Often for development, it is already enough to simply run the code and tests on
the host machine and/or rely on the CI for the tests. Sometimes, it is however
beneficial to run the full infrastructure on your development machine to
quickly iterate changes that rely on the full functionality, e.g. GUI, workers
or API. This section covers how to best setup the full infrasctructure for
development work.

Nomad depends on a set of databases, search engines, and other services. Those
must run to make use of nomad. We use *docker* and *docker-compose* to create a
unified environment that is easy to build and to run. During development we
may, however, wish to quickly reiterate the functionality of the API, gui or
workers. Because of this, these components are not run as a service within
docker during development work, and are instead run locally. This kind of setup
consists of the following steps:

1. Setup development dependencies on the host machine `as discussed above <#install-python-code-and-dependencies>`_
2. Run required infrastructure within docker containers that expose the services to the host machine:

    .. code-block:: sh

        cd nomad-FAIR/ops/docker-compose/infrastructure
        docker-compose -f docker-compose.yml -f docker-compose.prod.yml up -d mongo rabbitmq elastic

3. Run the API and worker on the host machine with the nomad cli tool that is installed in step 1:

    .. code-block:: sh

        nomad admin run appworker

4. Run the GUI on the host machine in development mode with yarn:

    .. code-block:: sh

        cd nomad-FAIR/gui
        yarn
        yarn start

Now you should have all the required services running. You can interact with
the application through the CLI, or through the GUI at localhost:3000. If you
started docker-compose in *deamon* mode (:code:`-d`) use :code:`docker-compose down` to
shut down the docker containers. Otherwise, just :code:`ctrl-c` the running output. 

.. note::
    The elastic service may require the host machine to increase the memory
    limit to function properly. `This is discussed here
    <https://www.elastic.co/guide/en/elasticsearch/reference/current/vm-max-map-count.html>`_.
    If you run into problems in starting the elastic service, you may need to
    increase the limit by running :code:`sysctl -w vm.max_map_count=262144` on the
    host machine.

Details on docker-compose for development
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You can use docker-compose to run all necessary databases with one single
docker-compose configuration. The :code:`docker-compose.yml` defines all the
different containers. We use docker-compose overrides to extend a base
configuration for different scenarios. The different overrides are:

- \*.prod.yml, production (to run the necessary databases for kubenetes deployments)
- \*.override.yml, development (development configuration, will be automatically used by docker-compose)
- \*.develk.yml, like development but also runs ELK

Running the tests
-----------------
You need to have the infrastructure partially running, as elaborated in :ref:`Build and run the development infrastructure with docker`.
The rest should be mocked or provided by the tests. Make sure that you do no
run any worker, as they will fight for tasks in the queue.

.. code-block:: sh

    cd nomad-FAIR/ops/docker-compose/infrastructure
    docker-compose -f docker-compose.yml -f docker-compose.prod.yml up -d rabbitmq elastic
    cd ../..
    pytest -svx tests

We use pylint, pycodestyle, and mypy to ensure code quality. To run those:

.. code-block:: sh

    nomad dev qa --skip-test

To run all tests and code qa:

.. code-block:: sh

    nomad dev qa

This mimiques the tests and checks that the GitLab CI/CD will perform.

Docker services for Nomad
-------------------------
As mentioned in the [introduction](setup.html#introduction), nomad comprises of
several interacting services. The services are used in various ways through
docker-compose files.  The purpose and name of each service stays the same
within different environments, but the details of the used image and
docker-compose configuration may vary across environments. For development,
only certain services are run, others are run directly by the host machine. For
production, everything runs within an a docker container. The following
subsections provide details on each service.

app
^^^
This service runs the nomad app, the API and the GUI. One can also run the
workers within this service, which can become useful in simple development
environments.

Before building the image, make sure to execute

.. code-block:: sh

    ./gitinfo.sh

This allows the app to present some information about the current git revision without
having to copy the git itself to the docker build context.

worker
^^^^^^
This service runs the celery worker that accepts jobs that are deployed by the
app service through the CLI or API.

mongo
^^^^^
Runs mongodb that is used for persisting worker states and other various
details about the repository and archive contents. On development environments
the mongodb ports are exposed to the host machine and you can access mongodb
via your preferred tools. Just make sure to use the right ports.

elastic
^^^^^^^
Runs ElasticSearch that provdes the search functionality. On development
environments the ElasticSearch ports are exposed to the host machine and you
can access mongodb via your preferred tools. Just make sure to use the right
ports.

rabbitmq
^^^^^^^^
Used as a message broker for celery tasks that are run by the worker service.

keycloak
^^^^^^^^
Runs keycloak for centralized user management.

Run nomad services
------------------
API and worker
^^^^^^^^^^^^^^
To simply run a worker with the installed nomad cli, do (from the root)

.. code-block:: sh

    nomad admin run worker

To run it directly with celery, do (from the root)

.. code-block:: sh

    celery -A nomad.processing worker -l info

You can also run worker and app together:

.. code-block:: sh

    nomad admin run appworker

Setup your (I)DE
----------------
The documentation section on development guidelines details how the code is organized,
tested, formatted, and documented. To help you meet these guidelines, we recomment to
use a proper IDE for development and ditch any VIM/Emacs (mal-)practices.

### Visual Studio Code

Here are some VSCode settings that will enable features for linting, some auto formating,
line size ruler, etc.

.. code-block:: json

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

Here are some example launch configs for VSCode:

.. code-block:: json

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
