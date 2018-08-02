This project tries and test approaches that might lead to an improved architecture for NOMAD (XT).

## Getting started

### Install the python in your own virtual environment.

```
pip install -r requirements.txt
pip install -e .
```

### Run dev infrastructure with docker.
```
cd ./infrastructure
docker-compose build
docker-compose up
```

Optionally register the infrastructue minio host to the minio client (mc).
```
mc config host add minio http://localhost:9007 AKIAIOSFODNN7EXAMPLE wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY
```

### Run the celery worker (should be moved to docker TODO)

```
celery -A nomad.processing worker -l info
```
You can use different debug level (e.g. switch `info` to `debug`)

Use watchdog during development. Install (i.e. [fixed](https://github.com/gorakhargosh/watchdog/issues/330) version fo MacOS)
```
pip install git+https://github.com/gorakhargosh/watchdog.git
```

Now use this to auto relead worker:
```
watchmedo auto-restart -d ./nomad -p '*.py' -- celery worker -l info -A nomad.processing
```

### Run tests.
```
python tests/test_files.py
```

### Generate the documentation

Install [pdoc](https://pypi.org/project/pdoc/) if necessary, and run this from the project root:
```
pdoc --html --html-dir=./build/docs --overwrite nomad
```

Do serve the current docs via HTTP during development, run
```
pdoc --http --http-port 8080
```
and open [http://localhost:8080/nomad](http://localhost:8080/nomad).

We will probably move to Sphynx at an apropriate moment.

## Contributing

### Code quality

- Use an IDE (e.g. [vscode](https://code.visualstudio.com/)) to enforce code [formatting and linting](https://code.visualstudio.com/docs/python/linting).
- There is a style guide to python. Write [pep-8](https://www.python.org/dev/peps/pep-0008/) compliant python code. An exception is the line cap at 79, which can be broken but keep it 90-ish.
- Use docstrings and document any *public* API of each submodule (e.g. python file). Public meaning API that is exposed to other submodules (i.e. other python files). Keep docstrings compliant to [pep-257](https://www.python.org/dev/peps/pep-0257/). Use sphynx references a lot!
- Be [pythonic](https://docs.python-guide.org/writing/style/) and watch [this](https://www.youtube.com/watch?v=wf-BqAjZb8M).
- The project structure is according to [this](https://docs.python-guide.org/writing/structure/) guide. Keep it!
- Test the public API of each submodule (i.e. python file)



