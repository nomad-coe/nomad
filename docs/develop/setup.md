---
title: Getting started
---
# Setup a dev environment, run the app, and run test

This is a step-by-step guide to get started with NOMAD development. You will clone
all sources, set-up a *Python* and *node* environment, install all necessary dependency,
run the infrastructure in development mode, learn to run out test-suites, and setup-up
*Visual Studio Code* for NOMAD development.

This is not about working with the NOMAD Python package. You can find the `nomad-lab`
documentation [here](../pythonlib.md).

## Clone the sources
If not already done, you should clone nomad. If you have a gitlab@MPCDF account, you can clone with git URL:

```
git clone git@gitlab.mpcdf.mpg.de:nomad-lab/nomad-FAIR.git nomad
```

Otherwise, clone using HTTPS URL:

```
git clone https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR.git nomad
```

then change directory to nomad

```
cd nomad
```

There are several branches in the repository. The master branch contains the latest released version,
but there is also a develop (new features) and release branch (hotfixes). There are also
tags for each version called vX.X.X. Checkout the branch you want to work on.
```
git checkout develop
```
The development branches are protected and you should create a new branch including your changes.
```
git checkout -b <my-branch-name>
```
This branch can be pushed to the repo, and then later may be merged to the relevant branch.

### Install sub-modules

Nomad is based on python modules from the NOMAD-coe project.
This includes parsers, python-common and the meta-info. These modules are maintained as
their own GITLab/git repositories. To clone and initialize them run:

```sh
git submodule update --init
```

## Installation

### Setup a Python environment

You should work in a Python virtual environment.

#### pyenv
The nomad code currently targets python 3.9. If your host machine has an older version installed,
you can use [pyenv](https://github.com/pyenv/pyenv) to use python 3.9 in parallel to your
system's python.

#### virtualenv
Create a virtual environment. It allows you
to keep nomad and its dependencies separate from your system's python installation.
Make sure that the virtual environment is based on Python 3.8 or higher (ideally Python 3.9). Use the built-in `venv` or (virtualenv)[https://pypi.org/project/virtualenv/] alternatively.

```
python3 -m venv .pyenv
source .pyenv/bin/activate
```

#### conda
If you are a conda user, there is an equivalent, but you have to install pip and the
right python version while creating the environment.
```sh
conda create --name nomad_env pip python=3.9
conda activate nomad_env
```

To install libmagick for conda, you can use (other channels might also work):
```sh
conda install -c conda-forge --name nomad_env libmagic
```

#### Upgrade pip
Make sure you have the most recent version of pip:
```sh
pip install --upgrade pip
```

### Install missing system libraries (e.g. on MacOS)

Even though the NOMAD infrastructure is written in python, there is a C library
required by one of our python dependencies. Libmagic is missing on some systems.
Libmagic allows to determine the MIME type of files. It should be installed on most
unix/linux systems. It can be installed on MacOS with homebrew:

```sh
brew install libmagic
```

If you are using an Mac with Apple Silicon, we recommend to use rosetta, homebrew
for Intel, and install and use an Intel based Python. The second answer in this
[Stackoverflow post](https://stackoverflow.com/questions/64882584/how-to-run-the-homebrew-installer-under-rosetta-2-on-m1-macbook)
describes how to use both the Apple and Intel homebrew simultaneously.

### Install nomad
The following command can be used to install all dependencies of all submodules and nomad
itself. If successful you can skip the rest of this *Install nomad* section.
```
./scripts/setup_dev_env.sh
```

Install all the requirements needed for development (including submodul requirements):

```sh
pip install --prefer-binary -r requirements-dev.txt
```

Finally, you can add nomad to the environment itself (including all extras).
The `-e` option will install the NOMAD with symbolic links allowing you
to change the code without having to reinstall after each change.
```sh
pip install -e .[parsing,infrastructure,dev]
```

If pip tries to use and compile sources and this creates errors, it can be told to prefer binary version:

```sh
pip install -e .[parsing,infrastructure,dev] --prefer-binary
```


### Update GUI artifacts

The NOMAD GUI requires static artifacts that are generated from the NOMAD Python codes.
```sh
python -m nomad.cli dev gui-artifacts --output-directory gui/src
python -m nomad.cli dev gui-config >gui/public/env.js
```

Or simply run
```sh
./scripts/generate_gui_artifacts.sh
```

The generated files are stored in GIT. The GUI code might not match the expected data in
outdated files. If there are changes to units, metainfo, new parsers, new toolkit notebooks it
might be necessary to regenerate these gui artifacts.

In addition, you have to do some more steps to prepare your working copy to run all
the tests. See below.

## Run the infrastructure

### Install docker
You need to install [docker](https://docs.docker.com/engine/install/).
Docker nowadays comes with `docker compose` build in. Prior, you needed to
install the standalone [docker-compose](https://docs.docker.com/compose/install/).

### Run required 3rd party services

To run NOMAD, some 3rd party services are needed

- elastic search: nomad's search and analytics engine
- mongodb: used to store processing state
- rabbitmq: a task queue used to distribute work in a cluster

All 3rd party services should be run via *docker-compose* (see below).
Keep in mind the *docker-compose* configures all services in a way that mirror
the configuration of the python code in `nomad/config.py` and the gui config in
`gui/.env.development`.

The default virtual memory for Elasticsearch is likely to be too low. On Linux, you can run the following command as root:
```sh
sysctl -w vm.max_map_count=262144
```

To set this value permanently, see [here](https://www.elastic.co/guide/en/elasticsearch/reference/current/vm-max-map-count.html). Then, you can run all services with:
```sh
cd ops/docker-compose/infrastructure
docker compose up -d mongo elastic rabbitmq
cd ../../..
```

If your system almost ran out of disk space the elasticsearch enforces a read-only index block ([read more](https://www.elastic.co/guide/en/elasticsearch/reference/6.2/disk-allocator.html)), but
after clearing up the disk space you need to reset it manually using the following command:

```sh
curl -XPUT -H "Content-Type: application/json" http://localhost:9200/_all/_settings -d '{"index.blocks.read_only_allow_delete": false}'
```

Note that the ElasticSearch service has a known problem in quickly hitting the
virtual memory limits of your OS. If you experience issues with the
ElasticSearch container not running correctly or crashing, try increasing the
virtual memory limits as shown [here](https://www.elastic.co/guide/en/elasticsearch/reference/current/vm-max-map-count.html).

To shut down everything, just `ctrl-c` the running output. If you started everything
in *deamon* mode (`-d`) use:
```sh
docker compose down
```

Usually these services are used only by NOMAD, but sometimes you also
need to check something or do some manual steps. You can access mongodb and elastic search
via your preferred tools. Just make sure to use the right ports.

## Run NOMAD

Before you run NOMAD for development purposes, you should configure it to use the `test`
realm of our user management system. By default, NOMAD will use the `fairdi_nomad_prod` realm.
Create a `nomad.yaml` file in the root folder:

```
keycloak:
  realm_name: fairdi_nomad_test
```

### App and Worker
NOMAD consist of the NOMAD app/api, a worker, and the GUI. You can run the app and the worker with
the NOMAD cli. These commands will run the services and display their log output. You should open
them in separate shells as they run continuously. They will not watch code changes and
you have to restart manually.

```sh
nomad admin run app
```

```sh
nomad admin run worker
```

Or both together in one process:
```
nomad admin run appworker
```

On MacOS you might run into multiprocessing errors. That can be solved as described [here](https://stackoverflow.com/questions/50168647/multiprocessing-causes-python-to-crash-and-gives-an-error-may-have-been-in-progr).

The app will run at port 8000 by default.

To run the worker directly with celery, do (from the root)
```sh
celery -A nomad.processing worker -l info
```

Before you can run the gui, make sure that generated artifacts are up-to-date:
```sh
./scripts/generate_gui_artifacts.sh
```

Also, make sure that the config file (`gui/public/env.js`) have been properly created:
```sh
./scripts/generate_gui_config.sh
```

If you run the gui on its own (e.g. with react dev server below), you also have to have
the app manually. The gui and its dependencies run on [node](https://nodejs.org) and
the [yarn](https://yarnpkg.com/) dependency manager. Read their documentation on how to
install them for your platform.
```sh
cd gui
yarn
yarn start
```

### JupyterHUB

NOMAD also has a build in JupyterHUB that is used to launch remote tools (e.g. Jupyter
notebooks).

To run the JupyterHUB, some additional configuration might be necessary.
```sh
north:
    hub_connect_ip: 'host.docker.internal'
    jupyterhub_crypt_key: '<crypt key>'
```

On Windows system, you might have to activate further specific functionality:
```sh
north:
    hub_connect_ip: 'host.docker.internal'
    hub_connect_url: 'http://host.docker.internal:8081'
    windows: true
    jupyterhub_crypt_key: '<crypt key>'
```

- If you are not on Linux, you need to configure how JupyterHUB can reach your host network from
docker containers. For Windows and MacOS you need to set `hub_connect_ip` to `host.docker.internal`. For linux you can leave it out and use the default `172.17.0.1`, unless you changed your
docker configuration.
- You have to generate a `crypt key` with `openssl rand -hex 32`.
- You might need to install [configurable-http-proxy](https://github.com/jupyterhub/configurable-http-proxy).

The *configurable-http-proxy* It comes as a node package. See [node](https://nodejs.org) for how to install `npm`. The proxy can be globally installed with:

```sh
npm install -g configurable-http-proxy
```

The JupyterHUB is a separate application. You can run the JuypterHUB similar
tp the other part.

```sh
nomad admin run hub
```

To run the JupyterHUB directly, do (from the root)
```sh
jupyterhub -f nomad/jupyterhub_config.py --port 9000
```

## Running tests

### Backend tests

To run the tests some additional settings and files are necessary that are not part
of the code base.

You have to provide static files to serve the docs and NOMAD distribution:
```sh
./scripts/generate_docs_artifacts.sh
rm -rf site && mkdocs build && mv site nomad/app/static/docs
```

You need to have the infrastructure partially running: elastic, rabbitmq.
The rest should be mocked or provided by the tests. Make sure that you do no run any
worker, as they will fight for tasks in the queue.
```sh
cd ops/docker-compose/infrastructure
docker compose up -d elastic rabbitmq
cd ../..
pytest -svx tests
```

We use pylint, pycodestyle, and mypy to ensure code quality. To run those:
```sh
nomad dev qa --skip-test
```

To run all tests and code qa:
```sh
nomad dev qa
```

This mimics the tests and checks that the GitLab CI/CD will perform.

### Frontend tests

We use
[`testing-library`](https://testing-library.com/docs/react-testing-library/intro/)
to implement our GUI tests and testing-library itself uses
[`jest`](https://jestjs.io/) to run the tests. Tests are written in `\*.spec.js`
files that accompany the implementation. Tests should focus on functionality,
not on implementation details: `testing-library` is designed to enforce this kind
of testing.

!!! note

    When testing HTML output, the elements are rendered using
    [jsdom](https://github.com/jsdom/jsdom): this is not completely identical
    to using an actual browser (does not support e.g. WebGL), but in practice
    is realistic enough for the majority of the test.

We have adopted a `pytest`-like structure for organizing the test utilities:
each source code folder may contain a `conftest.js` file that contains
utilities that are relevant for testing the code in that particular folder.
These utilities can usually be placed into the following categories:

 - Custom renders: When testing React components, the
   [`render`](https://testing-library.com/docs/react-testing-library/api/#render)-function
   is used to display them on the test DOM. Typically your components require
   some parts of the infrastructure to work properly, which is achieved by
   wrapping your component with other components that provide a context. Custom
   render functions can do this automatically for you. E.g. the default render
   as exported from `src/components/conftest.js` wraps your components with an
   infrastructure that is very similar to the production app. See
   [here](https://testing-library.com/docs/react-testing-library/setup/#custom-render)
   for more information.
 - Custom queries: See
   [here](https://testing-library.com/docs/react-testing-library/setup/#add-custom-queries)
   for more information.
 - Custom expects: These are reusable functions that perform actual tests using
   the expect-function. Whenever the same tests are performed by several
   \*.spec.js files, you should formalize these common tests into a
   `expect*`-function and place it in a relevant conftest.js file.

Often you components will need to communicate with the API during tests. One
should generally avoid using manually created mocks for the API traffic, and
instead prefer using API responses that originate from an actual API call
during testing. Manually created mocks require a lot of manual work in creating
them and keeping them up to date and true integration tests are impossible
to perform without live communication with an API. In order to simplify the API
communication during testing, you can use the `startAPI`+`closeAPI` functions, that
will prepare the API traffic for you. A simple example could look like this:

```javascript
import React from 'react'
import { waitFor } from '@testing-library/dom'
import { startAPI, closeAPI, screen } from '../../conftest'
import { renderSearchEntry, expectInputHeader } from '../conftest'

test('periodic table shows the elements retrieved through the API', async () => {
  startAPI('<state_name>', '<snapshot_name>')
  renderSearchEntry(...)
  expect(...)
  closeAPI()
})
```

Here the important parameters are:

 - `<state_name>`: Specifies an initial backend configuration for this test. These
   states are defined as python functions that are stored in
   nomad-FAIR/tests/states, example given below. These functions may e.g.
   prepare several uploads entries, datasets, etc. for the test.
 - `<snapshot_name>`: Specifies a filepath for reading/recording pre-recorded API
   traffic.

An example of a simple test state could look like this:

```python
from nomad import infrastructure
from nomad.utils import create_uuid
from nomad.utils.exampledata import ExampleData

def search():
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username="test")
    data = ExampleData(main_author=main_author)
    upload_id = create_uuid()
    data.create_upload(upload_id=upload_id, published=True, embargo_length=0)
    data.create_entry(
        upload_id=upload_id,
        entry_id=create_uuid(),
        mainfile="test_content/test_entry/mainfile.json",
        results={
            "material": {"elements": ["C", "H"]},
            "method": {},
            "properties": {}
        }
    )
    data.save()
```
When running in the `test-integration` or `test-record` mode (see below), this
function will be executed in order to prepare the application backend. The
`closeAPI` function will handle cleaning the test state between successive
`startAPI` calls: it will completely wipe out MongoDB, ElasticSearch and the
upload files.

!!! note

    The tests are using the configuration specified in `gui/tests/nomad.yaml`, that
    specifies a separate database/filesystem config in order to prevent interacting with
    any other instances of NOMAD.

In order to control how the API traffic is handled, there are three main ways
for running the test suite, as configured in `package.json`:

 - `yarn test [<filename>]`: Runs the tests parallelly in an 'offline'
   mode: `startAPI` will use pre-recorded API snapshot files that are found in
   gui/tests.
 - `yarn test-integration filename>]`: Runs the tests serially and `startAPI`
   will forward any API traffic to a live API that is running locally.
 - `yarn test-record [<filename>]`: Runs the tests serially and `startAPI` will
   forward traffic to a live API that is running locally, additionally
   recording the traffic to the specified snapshot file.

!!! note

    Before running against a live API (`yarn test-integration` and `yarn
    test-record`), you need to boot up the infrastructure and ensure that the
    nomad package is available with the correct test configuration:

    1. Have the docker infrastructure running: `docker compose up`

    2. Have the `nomad appworker` running with the config found in
       `gui/tests/nomad.yaml`. This can be achieved e.g. with the command: `export
       NOMAD_CONFIG=gui/tests/nomad.yaml; nomad admin run appworker`

    3. Activate the correct python virtual environment before running the tests
       with yarn (yarn will run the python functions that prepare the state).

## Setup your IDE

The documentation section for development guidelines (see below) details how the code is organized,
tested, formatted, and documented. To help you meet these guidelines, we recommend to
use a proper IDE for development and ditch any VIM/Emacs (mal-)practices.

We strongly recommend that all developers use *visual studio code*, or *vscode* for short,
(this is a completely different product than *visual studio*). It is available for free
for all major platforms [here](https://code.visualstudio.com/download).

You should launch and run vscode directly from the projects root directory. The source
code already contains settings for vscode in the `.vscode` directory. The settings
contain the same setup for stylechecks, linter, etc. that is also used in our CI/CD
pipelines. 
In order to ractually use the these features you have to make sure that they are enabled 
in your own User settings:
```
    "python.linting.pycodestyleEnabled": true,
    "python.linting.pylintEnabled": true,
    "python.linting.mypyEnabled": true,
    "python.testing.pytestEnabled": true,
```


The settings also include a few launch configuration for vscode's debugger. You can create
your own launch configs in `.vscode/launch.json` (also in .gitignore).

The settings expect that you have installed a python environment at `.pyenv` as
described in this tutorial (see above).

We also provide developers a vscode extension which is designed to support nomad schema language.
One can generate the extension using the following command after nomad installation

```sh
  nomad dev vscode-extension -o <path to output>
```

this command generate an up-to-date extension folder namely `nomad-vscode`. You can either copy
this folder into vscode extensions folder `~/.vscode/extensions/` or create an installable package as follows

```sh
  sudo npm install -g vsce # if vsce is not installed
  cd ./nomad-vscode
  vsce package
```

then install the extension by drag the file `nomad-0.0.x.vsix` and drop it into the extension panel of the vscode.
