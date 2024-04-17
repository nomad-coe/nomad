# How to get started in development

This is a step-by-step guide to get started with NOMAD development. You will clone
all sources, set up a *Python* and *Node.js* environment, install all necessary dependencies,
run the infrastructure in development mode, learn to run the test suites, and set up
*Visual Studio Code* for NOMAD development.

This is not about working with the NOMAD Python package `nomad-lab`. You can find its
documentation [here](../programmatic/pythonlib.md).

## Clone the sources

If not already done, you should clone NOMAD. If you have an account at the
[MPDCF Gitlab](https://gitlab.mpcdf.mpg.de/){:target="_blank"}, you can clone with the SSH URL:

```shell
git clone git@gitlab.mpcdf.mpg.de:nomad-lab/nomad-FAIR.git nomad
```

Otherwise, clone using the HTTPS URL:

```shell
git clone https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR.git nomad
```

Then change directory to `nomad`

```shell
cd nomad
```

There are several branches in the repository. The `master` branch contains the latest released version,
but there is also a `develop` (new features) and `release` branch (hotfixes). There are also
tags for each version called `vX.X.X`. Check out the branch you want to work on.

```shell
git checkout develop
```

The development branches are protected and you should create a new branch including your changes.

```shell
git checkout -b <my-branch-name>
```

This branch can be pushed to the repo, and then later may be merged to the relevant branch.

### Install submodules

Parts of the NOMAD software, such as parsers, are maintained in separate Git repositories.
These are then connected to the main repository as Git submodules. To clone and initialize
these submodules, run:

```shell
git submodule update --init
```

## Installation

### Set up a Python environment

The NOMAD code currently targets Python 3.9. You should work in a Python virtual environment.

#### Pyenv

If your host machine has an older version installed,
you can use [pyenv](https://github.com/pyenv/pyenv){:target="_blank"} to use Python 3.9 in parallel with your
system's Python.

#### Virtualenv

Create a virtual environment. It allows you
to keep NOMAD and its dependencies separate from your system's Python installation.
Make sure that the virtual environment is based on Python 3.9.
Use either the built-in `venv` module (see example) or [virtualenv](https://pypi.org/project/virtualenv/){:target="_blank"}.

```shell
python3 -m venv .pyenv
source .pyenv/bin/activate
```

#### Conda

If you are a conda user, there is an equivalent, but you have to install `pip` and the
right Python version while creating the environment.

```shell
conda create --name nomad_env pip python=3.9
conda activate nomad_env
```

To install libmagick for Conda, you can use (other channels might also work):

```shell
conda install -c conda-forge --name nomad_env libmagic
```

#### Upgrade pip

Make sure you have the most recent version of `pip`:

```shell
pip install --upgrade pip
```

### Install missing system libraries (e.g. on MacOS)

Even though the NOMAD infrastructure is written in Python, there is a C library
required by one of our Python dependencies. Libmagic is missing on some systems.
Libmagic allows to determine the MIME type of files. It should be installed on most
Unix/Linux systems. It can be installed on MacOS with homebrew:

```shell
brew install libmagic
```

If you are using a Mac with Apple Silicon, we recommend that you use rosetta, homebrew
for Intel, and install and use an Intel-based Python. The second answer in this
[Stackoverflow post](https://stackoverflow.com/questions/64882584/how-to-run-the-homebrew-installer-under-rosetta-2-on-m1-macbook){:target="_blank"}
describes how to use both the Apple and Intel homebrew simultaneously.

### Install NOMAD

The following command can be used to install all dependencies of all submodules
and NOMAD itself.

```shell
./scripts/setup_dev_env.sh
```

??? note "Installation details"
    Here is more detailed rundown of the installation steps.

    First we ensure that all submodules are up-to-date:

    ```shell
    git submodule update --init --recursive
    ```

    Previous build is cleaned:

    ```shell
    rm -rf nomad/app/static/docs
    rm -rf nomad/app/static/gui
    rm -rf site
    ```

    All the requirements needed for development (including submodule requirements) are installed:

    ```shell
    pip install --prefer-binary -r requirements-dev.txt
    ```

    Next we install the `nomad` package itself (including all extras). The `-e`
    option will install NOMAD with symbolic links that allow you to change
    the code without having to reinstall after each change.

    ```shell
    pip install -e .[parsing,infrastructure,dev]
    ```

    If pip tries to use and compile sources that create errors, it can be told to prefer the binary version:

    ```shell
    pip install -e .[parsing,infrastructure,dev] --prefer-binary
    ```

    The NOMAD GUI requires a static `.env` file, which can be generated with:

    ```shell
    python -m nomad.cli dev gui-env > gui/.env.development
    ```

    This file includes some of the server details needed so that the
    GUI can make the initial connection properly. If, for example, you change the server
    address in your NOMAD configuration file, it will be necessary to regenerate
    this `.env` file. In production this file will be overridden.

In addition, you have to do some more steps to prepare your working copy to run
all the tests, see below.

## Run the infrastructure

### Install Docker

You need to install [Docker](https://docs.docker.com/engine/install/){:target="_blank"}.
Docker nowadays comes with Docker Compose (`docker compose`) built-in. Prior, you needed to
install the standalone [Docker Compose (`docker-compose`)](https://docs.docker.com/compose/install/){:target="_blank"}.

### Run required 3rd party services

To run NOMAD, some 3rd party services are needed

- Elasticsearch: NOMAD's search and analytics engine
- MongoDB: used to store processing state
- RabbitMQ: a task queue used to distribute work in a cluster

All 3rd party services should be run via `docker compose` (see below).
Keep in mind that `docker compose` configures all services in a way that mirrors
the configuration of the Python code in `nomad/config.py` and the GUI config in
`gui/.env.development`.

The default virtual memory for Elasticsearch will likely be too low. On Linux, you can run the following command as root:

```shell
sysctl -w vm.max_map_count=262144
```

To set this value permanently, see [here](https://www.elastic.co/guide/en/elasticsearch/reference/current/vm-max-map-count.html){:target="_blank"}. Then you can run all services with:

```shell
cd ops/docker-compose/infrastructure
docker compose up -d elastic mongo rabbitmq
cd ../../..
```

If your system almost ran out of disk space, Elasticsearch enforces a read-only index block ([read more](https://www.elastic.co/guide/en/elasticsearch/reference/6.2/disk-allocator.html){:target="_blank"}), but
after clearing up the disk space you need to reset it manually using the following command:

```shell
curl -XPUT -H "Content-Type: application/json" http://localhost:9200/_all/_settings -d '{"index.blocks.read_only_allow_delete": false}'
```

Note that the Elasticsearch service has a known problem in quickly hitting the
virtual memory limits of your OS. If you experience issues with the
Elasticsearch container not running correctly or crashing, try increasing the
virtual memory limits as shown [here](https://www.elastic.co/guide/en/elasticsearch/reference/current/vm-max-map-count.html){:target="_blank"}.

To shut down everything, just `ctrl-c` the running output. If you started everything
in *deamon* mode (`-d`) use:

```shell
docker compose down
```

Usually these services are used only by NOMAD, but sometimes you also
need to check something or do some manual steps. You can access MongoDB and Elasticsearch
via your preferred tools. Just make sure to use the right ports.

## Run NOMAD

### `nomad.yaml`

Before you run NOMAD for development purposes, you should configure it to use the test
realm of our user management system. By default, NOMAD will use the `fairdi_nomad_prod` realm.
Create a `nomad.yaml` file in the root folder:

```yaml
keycloak:
  realm_name: fairdi_nomad_test
```

You might also want to exclude some of the default plugins, or only include the plugins
you'll need. Especially plugins with slower start-up and import times due to instantiation
of large schemas (e.g. nexus create couple thousand definitions for 70+ applications) can
often be excluded:

```yaml
plugins:
  exclude:
    - parsers/nexus
```

Note that this will lead to [failing tests](#backend-tests) for the excluded plugins.

### App and worker

NOMAD consists of the NOMAD app/API, a worker, and the GUI. You can run the app and the worker with
the NOMAD CLI. These commands will run the services and display their log output. You should open
them in separate shells as they run continuously. They will not watch code changes and
you have to restart manually.

```shell
nomad admin run app
```

```shell
nomad admin run worker
```

Or both together in one process:

```shell
nomad admin run appworker
```

On MacOS you might run into multiprocessing errors. That can be solved as described [here](https://stackoverflow.com/questions/50168647/multiprocessing-causes-python-to-crash-and-gives-an-error-may-have-been-in-progr){:target="_blank"}.

The app will run at port 8000 by default.

To run the worker directly with Celery, do (from the root)

```shell
celery -A nomad.processing worker -l info
```

If you run the GUI on its own (e.g. with the React dev server below), you also need to start
the app manually. The GUI and its dependencies run on [Node.js](https://nodejs.org){:target="_blank"} and
the [Yarn](https://yarnpkg.com/){:target="_blank"} dependency manager. Read their documentation on how to
install them for your platform.

```shell
cd gui
yarn
yarn start
```

Note that the current codebase requires Node.js version 16.20. If you have a newer version installed on your system, you may need to downgrade or use the following workaround as described [here](https://codingbeautydev.com/blog/node-err-ossl-evp-unsupported/) before `yarn start` works properly.
```
export NODE_OPTIONS=--openssl-legacy-provider
```

### JupyterHub

NOMAD also has a built-in JupyterHub that is used to launch remote tools (e.g. Jupyter
notebooks).

To run JupyterHub, some additional configuration might be necessary.

```yaml
north:
  hub_connect_ip: 'host.docker.internal'
  jupyterhub_crypt_key: '<crypt key>'
```

On Windows system, you might have to activate further specific functionality:

```yaml
north:
  hub_connect_ip: 'host.docker.internal'
  hub_connect_url: 'http://host.docker.internal:8081'
  windows: true
  jupyterhub_crypt_key: '<crypt key>'
```

- If you are not on Linux, you need to configure how JupyterHub can reach your host
  network from docker containers. For Windows and MacOS you need to set `hub_connect_ip`
  to `host.docker.internal`. For Linux you can leave it out and use the default
  `172.17.0.1`, unless you changed your docker configuration.

- You have to generate a `crypt key` with `openssl rand -hex 32`.

- You might need to install
  [configurable-http-proxy](https://github.com/jupyterhub/configurable-http-proxy){:target="_blank"}.

The `configurable-http-proxy` comes as a Node.js package. See
[Node.js](https://nodejs.org){:target="_blank"} for how to install `npm`. The proxy can be globally
installed with:

```shell
npm install -g configurable-http-proxy
```

JupyterHub is a separate application. You can run JuypterHub similar to the other part:

```shell
nomad admin run hub
```

To run JupyterHub directly, do (from the root)

```shell
jupyterhub -f nomad/jupyterhub_config.py --port 9000
```

## Running tests

### Backend tests

To run the tests some additional settings and files are necessary that are not part
of the codebase.

You have to provide static files to serve the docs and NOMAD distribution:

```shell
./scripts/generate_docs_artifacts.sh
rm -rf site && mkdocs build && mv site nomad/app/static/docs
```

You need to have the infrastructure partially running: `elastic`, `mongo`, `rabbitmq`.
The rest should be mocked or provided by the tests. Make sure that you do not run any
worker, as they will fight for tasks in the queue. To start the infrastructure and run the
tests, use:

```shell
cd ops/docker-compose/infrastructure
docker compose up -d elastic mongo rabbitmq
cd ../../..
pytest -sv tests
```

!!! note
    Some of these tests will fail because a few large files are not included in the Git
    repository. You may ignore these for local testing, they are still checked by the
    CI/CD pipeline:

    ```text
    FAILED tests/archive/test_archive.py::test_read_springer - AttributeError: 'NoneType' object has no attribute 'seek'
    FAILED tests/normalizing/test_material.py::test_material_bulk - assert None
    FAILED tests/normalizing/test_system.py::test_springer_normalizer - IndexError: list index out of range
    ```

    If you excluded plugins in your [NOMAD config](### `nomad.yaml`), then those tests
    will also fail.

We use Ruff and Mypy to maintain code quality. Additionally, we recommend installing the Ruff [plugins](https://docs.astral.sh/ruff/integrations/){:target="_blank"} for your code editor to streamline the process. To execute Ruff and Mypy from the command line, you can utilize the following command:
```shell
nomad dev qa --skip-tests
```

We use ruff as a linter and as an autoformatter. If you only want to lint your code, you can run:
```shell
ruff check .
```

To format your code you can run:
```shell
ruff format .
```

To run all tests and code QA:

```shell
nomad dev qa
```

This mimics the tests and checks that the GitLab CI/CD will perform.

If you are migrating an old merge request to a formatted one, please see find the migration guide [here](./migrate-to-autoformatter.md)

#### Custom pytest options

##### `--celery-inspect-timeout`

To ensure that all tests are independent despite reusing the same queue and workers, the `worker` fixture cleans up all running tasks after the test. This involves by default a timeout of one second, which accumulates over all tests using that fixture. In one local development environment this made a difference of about 14 min vs. 7 min without the timeout.

If you want to speed up your local testing, you can use `pytest --celery-inspect-timeout 0.1` to shorten the timeout to a tenth of a second. Be aware that this might leave tasks running, which can affect later tests.

##### `--fixture-filters`

You may want to run only tests that use a specific fixture, e.g. if you are editing that fixture. For example, to run only tests with the `worker` fixture, use `pytest --fixture-filters worker`. If you list more than one, all of them must be requested by the test to be included.

You can also negate a fixture by prefixing its name with `!`. For example, to run all tests that do not depend on the `worker` fixture, use `pytest --fixture-filters '!worker'` (quotes are needed for `!`).

Note that if `test1` depends on `fixture1`, and `fixture1` depends on `fixture2`, then `test1` also depends on `fixture2`, even though it was not explicitly listed as a parameter.

### Frontend tests

We use
[`testing-library`](https://testing-library.com/docs/react-testing-library/intro/){:target="_blank"}
to implement our GUI tests and `testing-library` itself uses
[`Jest`](https://jestjs.io/){:target="_blank"} to run the tests. Tests are written in `*.spec.js`
files that accompany the implementation. Tests should focus on functionality,
not on implementation details: `testing-library` is designed to enforce this kind
of testing.

!!! note

    When testing HTML output, the elements are rendered using
    [jsdom](https://github.com/jsdom/jsdom){:target="_blank"}: this is not completely identical
    to using an actual browser (e.g. does not support WebGL), but in practice
    is realistic enough for the majority of the test.

#### Test structure

We have adopted a `pytest`-like structure for organizing the test utilities:
each source code folder may contain a `conftest.js` file that contains
utilities that are relevant for testing the code in that particular folder.
These utilities can usually be placed into the following categories:

- Custom renders: When testing React components, the
  [`render`](https://testing-library.com/docs/react-testing-library/api/#render){:target="_blank"} function
  is used to display them on the test DOM. Typically your components require
  some parts of the infrastructure to work properly, which is achieved by
  wrapping your component with other components that provide a context. Custom
  render functions can do this automatically for you, e.g. the default render
  as exported from `src/components/conftest.js` wraps your components with an
  infrastructure that is very similar to the production app. See
  [here](https://testing-library.com/docs/react-testing-library/setup/#custom-render){:target="_blank"}
  for more information.

- Custom queries: See
  [here](https://testing-library.com/docs/react-testing-library/setup/#add-custom-queries){:target="_blank"}
  for more information.

- Custom expects: These are reusable functions that perform actual tests using
  the `expect` function. Whenever the same tests are performed by several
  `*.spec.js` files, you should formalize these common tests into an
  `expect*` function and place it in a relevant `conftest.js` file.

Often your components will need to communicate with the API during tests. One
should generally avoid using manually created mocks for the API traffic, and
instead prefer using API responses that originate from an actual API call
during testing. Manually created mocks require a lot of manual work in creating
them and keeping them up-to-date and true integration tests are impossible
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
  states are defined as Python functions that are stored in
  `nomad-FAIR/tests/states`, example given below. These functions may, for example,
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

When running in the online mode (see below), this function will be executed in
order to prepare the application backend. The `closeAPI` function will handle
cleaning the test state between successive `startAPI` calls: it will completely
wipe out MongoDB, Elasticsearch and the upload files.

#### Running frontend tests

The tests can be run in two different modes. *Offline testing* uses
pre-recorded files to mock the API traffic during testing. This allows one to
run tests more quickly without a server. During *online testing*, the tests
perform calls to a running server where a test state has been prepared. This
mode can be used to perform integration tests but also to record the snapshot
files needed by the offline testing.

##### Offline testing

This is the way our CI pipeline runs the tests and should be used locally, e.g.
whenever you wish to reproduce pipeline errors or when your tests do not
involve any API traffic.

1. Ensure that the GUI artifacts are up-to-date:

   ```shell
   ./scripts/generate_gui_test_artifacts.sh
   ```

   As snapshot tests do not connect to the server, the artifacts cannot be
   fetched dynamically from the server and static files need to be used
   instead.

2. Run `yarn test` to run the whole suite or `yarn test [<filename>]` to run a
   specific test.

##### Online testing

When you wish to record API traffic for offline testing, or to perform
integration tests, you will need to have a server running with the correct
configuration. To do this, follow these steps:

1. Have the docker infrastructure running: `docker compose up`

2. Have the `nomad appworker` running with the config found in
   `gui/tests/nomad.yaml`:
   `export NOMAD_CONFIG=gui/tests/nomad.yaml; nomad admin run appworker`

3. Activate the correct Python virtual environment before running the tests
   with Yarn (Yarn will run the Python functions that prepare the state).

4. Run the tests with `yarn test-record [<filename>]` if you wish to record a
   snapshot file or `yarn test-integration [<filename>]` if you want the
   perform the test without any recording.

## Build the Docker image

Normally the Docker image is build via a CI/CD pipeline that is run when pushing commits
to [NOMAD's GitLab at MPCDF](https://gitlab.mpcdf.mpg.de/){:target="_blank"}. These images are distributed
via NOMAD's GitLab container registry. For most purposes you would use these
automatically-built images.

If you want to build a custom image, e.g. to be used in your NOMAD Oasis, you can
run the NOMAD Docker build manually. From the cloned project root run:

```shell
docker build -t <image-name>:<image-tag> .
```

This will build the normal image intended for production use. There are other build
targets: `dev_python` and `dev_node`. Especially `dev_python` might be interesting
for debugging purposes as it contains all sources and dev dependencies. You
can build specific targets with:

```shell
docker build --target dev_python -t <image-name>:<image-tag> .
```

If you want to build an image directly from a remote Git repository (e.g. for a specific `branch`), run:

```shell
DOCKER_BUILDKIT=1 docker build --build-arg BUILDKIT_CONTEXT_KEEP_GIT_DIR=1 --pull -t <image-name>:<image-tag> https://github.com/nomad-coe/nomad.git#<branch>
```

The BuildKit parametrization ensures that the `.git` directory is available in the
Docker build context. NOMAD's build process requires the `.git` folder to determine the
package version from version tags in the repository.

The build process installs a substantial amount of dependencies and requires multiple
Docker images for various build stages. Make sure that Docker has at least 20 GB of
storage available.

## Setup your IDE

The documentation section for development guidelines (see below) provide details on how
the code is organized, tested, formatted, and documented. To help you meet these
guidelines, we recommend to use a proper IDE for development and ditch any Vim/Emacs
(mal-)practices.

We strongly recommend that all developers use *Visual Studio Code (VS Code)*. (This is a
completely different product than *Visual Studio*.) It is available for free
for all major platforms [here](https://code.visualstudio.com/download){:target="_blank"}.

You should launch and run VS Code directly from the project's root directory. The source
code already contains settings for VS Code in the `.vscode` directory. The settings
contain the same setup for style checks, linter, etc. that is also used in our CI/CD
pipelines. In order to actually use the these features, you have to make sure that they
are enabled in your own User settings:

```json
"python.linting.mypyEnabled": true,
"python.testing.pytestEnabled": true,
"[python]": {
  "editor.formatOnSave": true,
  "editor.defaultFormatter": "charliermarsh.ruff"
}
```

The settings also include a few launch configuration for VS Code's debugger. You can create
your own launch configs in `.vscode/launch.json` (also in `.gitignore`).

The settings expect that you have installed a Python environment at `.pyenv` as
described in this tutorial (see above).
