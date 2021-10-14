# Developing NOMAD

## Getting started

### Clone the sources
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

There are several branches in the repository. The master branch contains the latest released version, but there are also
develop branches for each version called vX.X.X. Checkout the branch you want to work on it
```
git checkout vX.X.X
```
The development branches are protected and you should create a new branch including your changes.
```
git checkout -b <my-branch-name>
```
This branch can be pushed to the repo, and then later may be merged to the relevant branch.

### Prepare your Python environment

You work in a Python virtual environment.

#### pyenv
The nomad code currently targets python 3.7. If you host machine has an older version installed,
you can use [pyenv](https://github.com/pyenv/pyenv) to use python 3.7 in parallel to your
system's python. Never the less, we have good experience with 3.8 and 3.9 users as well
and everything might work with newer versions as well.

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

If you use *pyenv* (or similar solutions) make sure that the `-p` arguments evaluates
to the `python` binary with the desired version.

#### conda
If you are a conda user, there is an equivalent, but you have to install pip and the
right python version while creating the environment.
```sh
conda create --name nomad_env pip python=3.7
conda activate nomad_env
```

To install libmagick for conda, you can use (other channels might also work):
```sh
conda install -c conda-forge --name nomad_env libmagic
```

## Setup
Using the following command one can install all the dependencies, and the sub-modules from the NOMAD-coe project
```
bash setup.sh
```

The script includes the following steps:

### 1. pip
Make sure you have the most recent version of pip:
```sh
pip install --upgrade pip
```


#### Missing system libraries (e.g. on MacOS)

Even though the NOMAD infrastructure is written in python, there is a C library
required by one of our python dependencies. Libmagic is missing on some systems.
Libmagic allows to determine the MIME type of files. It should be installed on most
unix/linux systems. It can be installed on MacOS with homebrew:

```sh
brew install libmagic
```

### 2. Install sub-modules
Nomad is based on python modules from the NOMAD-coe project.
This includes parsers, python-common and the meta-info. These modules are maintained as
their own GITLab/git repositories. To clone and initialize them run:

```sh
git submodule update --init
```

All requirements for these submodules need to be installed and they need to be installed
themselves as python modules. Run the `dependencies.sh` script that will install
everything into your virtual environment:
```sh
./dependencies.sh -e
```

If one of the Python packages that are installed during this process, fails because it
cannot be compiled on your platform, you can try `pip install --prefer-binary <packagename>`
to install set package manually.

The `-e` option will install the NOMAD-coe dependencies with symbolic links allowing you
to change the downloaded dependency code without having to reinstall after.

### 3. Install nomad
Finally, you can add nomad to the environment itself (including all extras)
```sh
pip install -e .[all]
```

If pip tries to use and compile sources and this creates errors, it can be told to prefer binary version:

```sh
pip install -e .[all] --prefer-binary
```

### 4. Generate GUI artifacts
The NOMAD GUI requires static artifacts that are generated from the NOMAD Python codes.
```sh
nomad dev metainfo > gui/src/metainfo.json
nomad dev search-quantities > gui/src/searchQuantities.json
nomad dev units > gui/src/units.js
./gitinfo.sh
```

In additional, you have to do some more steps to prepare your working copy to run all
the tests. See below.

## Install docker
One needs to install [docker](https://docs.docker.com/get-docker/) and [docker-compose](https://docs.docker.com/compose/install/).

## Running the infrastructure

To run NOMAD, some 3-rd party services are needed
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
docker-compose up -d mongo elastic rabbitmq
cd ../../..
```

If your system almost ran out of disk space the elasticsearch enforces a read-only index block ([read more](https://www.elastic.co/guide/en/elasticsearch/reference/6.2/disk-allocator.html)), but
after clearing up the disk space you need to reset it manually using the following command:

```sh
curl -XPUT -H "Content-Type: application/json" http://localhost:9200/_all/_settings -d '{"index.blocks.read_only_allow_delete": false}'
```

To shut down everything, just `ctrl-c` the running output. If you started everything
in *deamon* mode (`-d`) use:
```sh
docker-compose down
```

Usually these services only used by NOMAD, but sometimes you also
need to check something or do some manual steps. You can access mongodb and elastic search
via your preferred tools. Just make sure to use the right ports.

## Running NOMAD

Before you run NOMAD for development purposes, you should configure it to use the `test`
realm of our user management system. By default, NOMAD will use the `fairdi_nomad_prod` realm.
Create a `nomad.yaml` file in the root folder:

```
keycloak:
  realm_name: fairdi_nomad_test
```

NOMAD consist of the NOMAD app/api, a worker, and the GUI. You can run app and worker with
the NOMAD cli. These commands will run the services and show their logout put. You should open
them in separate shells as they run continuously. They will not watch code changes and
you have to restart manually.

```sh
nomad admin run app
```

```sh
nomad admin run worker
```

Or both together in once process:
```
nomad admin run appworker
```

The app will run at port 8000 by default.

To run the worker directly with celery, do (from the root)
```sh
celery -A nomad.processing worker -l info
```

When you run the gui on its own (e.g. with react dev server below), you have to have
the app manually also. The gui and its dependencies run on [node](https://nodejs.org) and
the [yarn](https://yarnpkg.com/) dependency manager. Read their documentation on how to
install them for your platform.
```sh
cd gui
yarn
yarn start
```

## Running tests

To run the tests some additional settings and files are necessary that are not part
of the code base.

First, you need to provide the `springer.msg` Springer materials database. It can
be copied from `/nomad/fairdi/db/data/springer.msg` on our servers and should
be placed at `nomad/normalizing/data/springer.msg`.

Second, you have to provide static files to serve the docs and NOMAD distribution:
```sh
cd docs
make html
cd ..
python setup.py compile
python setup.py sdist
cp dist/nomad-lab-*.tar.gz dist/nomad-lab.tar.gz
```

You need to have the infrastructure partially running: elastic, rabbitmq.
The rest should be mocked or provided by the tests. Make sure that you do no run any
worker, as they will fight for tasks in the queue.
```sh
cd ops/docker-compose/infrastructure
docker-compose up -d elastic rabbitmq
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
    "editor.rulers": [90],
    "editor.renderWhitespace": "all",
    "editor.tabSize": 4,
    "[javascript]": {
        "editor.tabSize": 2
    },
    "files.trimTrailingWhitespace": true,
    "editor.codeActionsOnSave": ["source.fixAll.eslint"],
    "python.linting.pylintEnabled": true,
    "python.linting.pylintArgs": [
        "--load-plugins=pylint_mongoengine,nomad/metainfo/pylint_plugin",
    ],
    "python.linting.pycodestylePath": "pycodestyle",
    "python.linting.pycodestyleEnabled": true,
    "python.linting.pycodestyleArgs": ["--ignore=E501,E701,E731"],
    "python.linting.mypyEnabled": true,
    "python.linting.mypyArgs": [
        "--ignore-missing-imports",
        "--follow-imports=silent",
        "--no-strict-optional"
    ],
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

## Code guidelines

### Design principles

- simple first, complicated only when necessary
- adopting generic established 3rd party solutions before implementing specific solutions
- only uni directional dependencies between components/modules, no circles
- only one language: Python (except, GUI of course)

### Rules

The are some *rules* or better strong *guidelines* for writing code. The following
applies to all python code (and were applicable, also to JS and other code):

- Use an IDE (e.g. [vscode](https://code.visualstudio.com/) or otherwise automatically
  enforce code ([formatting and linting](https://code.visualstudio.com/docs/python/linting)).
  Use `nomad qa` before committing. This will run all tests, static type checks, linting, etc.

- There is a style guide to python. Write [pep-8](https://www.python.org/dev/peps/pep-0008/)
  compliant python code. An exception is the line cap at 79, which can be broken but keep it 90-ish.

- Test the public API of each sub-module (i.e. python file)

- Be [pythonic](https://docs.python-guide.org/writing/style/) and watch
  [this](https://www.youtube.com/watch?v=wf-BqAjZb8M).

- Document any *public* API of each sub-module (e.g. python file). Public meaning API that
  is exposed to other sub-modules (i.e. other python files).

- Use google [docstrings](http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html).

- Add your doc-strings to the sphinx documentation in `docs`. Use .md, follow the example.
  Markdown in sphinx is supported via [recommonmark]
  (https://recommonmark.readthedocs.io/en/latest/index.html#autostructify)
  and [AutoStructify](http://recommonmark.readthedocs.io/en/latest/auto_structify.html)

- The project structure is according to [this guide](https://docs.python-guide.org/writing/structure/).
  Keep it!

- Write tests for all contributions.


### Enforcing Rules: CI/CD


These *guidelines* are partially enforced by CI/CD. As part of CI all tests are run on all
branches; further we run a *linter*, *pep8* checker, and *mypy* (static type checker). You can
run `nomad qa` to run all these tests and checks before committing.

The CI/CD will run on all refs that do not start with `dev-`. The CI/CD will
not release or deploy anything automatically, but it can be manually triggered after the
build and test stage completed successfully.

### Names and identifiers

There are is some terminology consistently used in this documentation and the source
code. Use this terminology for identifiers.

Do not use abbreviations. There are (few) exceptions: `proc` (processing); `exc`, `e` (exception);
`calc` (calculation), `repo` (repository), `utils` (utilities), and `aux` (auxiliary).
Other exceptions are `f` for file-like streams and `i` for index running variables.
Btw., the latter is almost never necessary in python.

Terms:

- upload: A logical unit that comprises one (.zip) file uploaded by a user.
- calculation: A computation in the sense that is was created by an individual run of a CMS code.
- raw file: User uploaded files (e.g. part of the uploaded .zip), usually code input or output.
- upload file/uploaded file: The actual (.zip) file a user uploaded
- mainfile: The mainfile output file of a CMS code run.
- aux file: Additional files the user uploaded within an upload.
- repo entry: Some quantities of a calculation that are used to represent that calculation in the repository.
- archive data: The normalized data of one calculation in nomad's meta-info-based format.

Throughout nomad, we use different ids. If something
is called *id*, it is usually a random uuid and has no semantic connection to the entity
it identifies. If something is called a *hash* than it is a hash build based on the
entity it identifies. This means either the whole thing or just some properties of
said entities.

- The most common hashes is the `calc_hash` based on mainfile and auxfile contents.
- The `upload_id` is a UUID assigned at upload time and never changed afterwards.
- The `mainfile` is a path within an upload that points to a main code output file.
  Since, the upload directory structure does not change, this uniquely ids a calc within the upload.
- The `calc_id` (internal calculation id) is a hash over the `mainfile` and respective
  `upload_id`. Therefore, each `calc_id` ids a calc on its own.
- We often use pairs of `upload_id/calc_id`, which in many context allow to resolve a calc
  related file on the filesystem without having to ask a database about it.
- The `pid` or (`coe_calc_id`) is an sequential interger id.
- Calculation `handle` or `handle_id` are created based on those `pid`.
  To create hashes we use :py:func:`nomad.utils.hash`.


### Logging

There are three important prerequisites to understand about nomad-FAIRDI's logging:

- All log entries are recorded in a central elastic search database. To make this database
  useful, log entries must be sensible in size, frequence, meaning, level, and logger name.
  Therefore, we need to follow some rules when it comes to logging.
- We use an *structured* logging approach. Instead of encoding all kinds of information
  in log messages, we use key-value pairs that provide context to a log *event*. In the
  end all entries are stored as JSON dictionaries with `@timestamp`, `level`,
  `logger_name`, `event` plus custom context data. Keep events very short, most
  information goes into the context.
- We use logging to inform about the state of nomad-FAIRDI, not about user
  behavior, input, data. Do not confuse this when determining the log-level for an event.
  For example, a user providing an invalid upload file, for example, should never be an error.

Please follow the following rules when logging:

- If a logger is not already provided, only use
  :py:func:`nomad.utils.get_logger` to acquire a new logger. Never use the
  build-in logging directly. These logger work like the system loggers, but
  allow you to pass keyword arguments with additional context data. See also
  the [structlog docs](https://structlog.readthedocs.io/en/stable/).
- In many context, a logger is already provided (e.g. api, processing, parser, normalizer).
  This provided logger has already context information bounded. So it is important to
  use those instead of acquiring your own loggers. Have a look for methods called
  `get_logger` or attributes called `logger`.
- Keep events (what usually is called *message*) very short. Examples are: *file uploaded*,
  *extraction failed*, etc.
- Structure the keys for context information. When you analyse logs in ELK, you will
  see that the set of all keys over all log entries can be quit large. Structure your
  keys to make navigation easier. Use keys like `nomad.proc.parser_version` instead of
  `parser_version`. Use module names as prefixes.
- Don't log everything. Try to anticipate, how you would use the logs in case of bugs,
  error scenarios, etc.
- Don't log sensitive data.
- Think before logging data (especially dicts, list, numpy arrays, etc.).
- Logs should not be abused as a *printf*-style debugging tool.

The following keys are used in the final logs that are piped to Logstash.
Notice that the key name is automatically formed by a separate formatter and
may differ from the one used in the actual log call.

Keys that are autogenerated for all logs:

 - `@timestamp`: Timestamp for the log
 - `@version`: Version of the logger
 - `host`: The host name from which the log originated
 - `path`: Path of the module from which the log was created
 - `tags`: Tags for this log
 - `type`: The *message_type* as set in the LogstashFormatter
 - `level`: The log level: `DEBUG`, `INFO`, `WARNING`, `ERROR`
 - `logger_name`: Name of the logger
 - `nomad.service`: The service name as configured in `config.py`
 - `nomad.release`: The release name as configured in `config.py`

Keys that are present for events related to processing an entry:

 - `nomad.upload_id`: The id of the currently processed upload
 - `nomad.calc_id`: The id of the currently processed entry
 - `nomad.mainfile`: The mainfile of the currently processed entry

Keys that are present for events related to exceptions:

 - `exc_info`: Stores the full python exception that was encountered. All
   uncaught exceptions will be stored automatically here.
 - `digest`: If an exception was raised, the last 256 characters of the message
   are stored automatically into this key. If you wish to search for exceptions
   in Kibana, you will want to use this value as it will be indexed unlike the
   full exception object.


### Copyright Notices

We follow this [recommendation](https://www.linuxfoundation.org/blog/2020/01/copyright-notices-in-open-source-software-projects/)
of the Linux Foundation for the copyright notice that is placed on top of each source
code file.

It is intended to provide a broad generic statement that allows all authors/contributors
of the NOMAD project to claim their copyright, independent of their organization or
individual ownership.

You can simply copy the notice from another file. From time to time we can use a tool
like [licenseheaders](https://pypi.org/project/licenseheaders/) to ensure correct
notices. In addition we keep an purely informative AUTHORS file.


## Git/GitLab

### Branches and clean version history

The `master` branch of our repository is *protected*. You must not (even if you have
the rights) commit to it directly. The `master` branch references the latest official
release (i.e. what the current NOMAD runs on). The current development is represented by
*version* branches, named `vx.x.x`. Usually there are two or more of these branched,
representing the development on *minor/bugfix* versions and the next *major* version(s).
Ideally these *version* branches are also not manually push to.

Instead you develop
on *feature* branches. These are branches that are dedicated to implement a single feature.
They are short lived and only exist to implement a single feature.

The lifecycle of a *feature* branch should look like this:

- create the *feature* branch from the last commit on the respective *version* branch that passes CI

- do your work and push until you are satisfied and the CI passes

- create a merge request on GitLab

- discuss the merge request on GitLab

- continue to work (with the open merge request) until all issues from the discussion are resolved

- the maintainer performs the merge and the *feature* branch gets deleted

While working on a feature, there are certain practices that will help us to create
a clean history with coherent commits, where each commit stands on its own.

```sh
  git commit --amend
```

If you committed something to your own feature branch and then realize by CI that you have
some tiny error in it that you need to fix, try to amend this fix to the last commit.
This will avoid unnecessary tiny commits and foster more coherent single commits. With *amend*
you are basically adding changes to the last commit, i.e. editing the last commit. If
you push, you need to force it `git push origin feature-branch --force-with-lease`. So be careful, and
only use this on your own branches.

```sh
  git rebase <version-branch>
```

Lets assume you work on a bigger feature that takes more time. You might want to merge
the version branch into your feature branch from time to time to get the recent changes.
In these cases, use rebase and not merge. Rebase puts your branch commits in front of the
merged commits instead of creating a new commit with two ancestors. It basically moves the
point where you initially branched away from the version branch to the current position in
the version branch. This will avoid merges, merge commits, and generally leave us with a
more consistent history.  You can also rebase before create a merge request, basically
allowing for no-op merges. Ideally the only real merges that we ever have, are between
version branches.

```sh
  git merge --squash <other-branch>
```

When you need multiple branches to implement a feature and merge between them, try to
use *squash*. Squashing basically puts all commits of the merged branch into a single commit.
It basically allows you to have many commits and then squash them into one. This is useful
if these commits where just made for synchronization between workstations or due to
unexpected errors in CI/CD, you needed a save point, etc. Again the goal is to have
coherent commits, where each commits makes sense on its own.

Often a feature is also represented by an *issue* on GitLab. Please mention the respective
issues in your commits by adding the issue id at the end of the commit message: *My message. #123*.

We tag releases with `vX.X.X` according to the regular semantic versioning practices.
After releasing and tagging the *version* branch is removed. Do not confuse tags with *version* branches.
Remember that tags and branches are both Git references and you can accidentally pull/push/checkout a tag.

The main NOMAD GitLab-project (`nomad-fair`) uses Git-submodules to maintain its
parsers and other dependencies. All these submodules are places in the `/dependencies`
directory. There are helper scripts to install (`./dependencies.sh` and
commit changes to all submodules (`./dependencies-git.sh`). After merging or checking out,
you have to make sure that the modules are updated to not accidentally commit old
submodule commits again. Usually you do the following to check if you really have a
clean working directory.

```sh
  git checkout something-with-changes
  git submodule update
  git status
```

### Submodules

We currently use git submodules to manage NOMAD internal dependencies (e.g. parsers).
All dependencies are python packages and installed via pip to your python environement.

This allows us to target (e.g. install) individual commits. More importantly, we can address c
ommit hashes to identify exact parser/normalizer versions. On the downside, common functions
for all dependencies (e.g. the python-common package, or nomad_meta_info) cannot be part
of the nomad-FAIRDI project. In general, it is hard to simultaneously develop nomad-FAIRDI
and NOMAD-coe dependencies.

Another approach is to integrate the NOMAD-coe sources with nomad-FAIRDI. The lacking
availability of individual commit hashes, could be replaces with hashes of source-code
files.

We use the `master` branch on all dependencies. Of course feature branches can be used on
dependencies to manage work in progress.