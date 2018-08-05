# Setup

### Install the legacy NOMAD submoduels.
This has to be done differently in the future. For no init the submodules and checkout
working branches/tags:
- submodules/parsers/parser-vasp master
- submodules/python-common master
- submodules/nomad-meta-info 1.6.0

To checkout a tag use:
```
git fetch --all --tags --prune
git checkout tags/1.6.0 -b 1.6.0
```

`pip install -r requirements` in `python-common`, and `pip install -e .` in `python-common` and
`parsers/parser-vasp`. Futhermore, there are some dependency issues in `python-commons` requirments.

### Install the python in your own virtual environment.

```
pip install -r requirements.txt
pip install -e .
```

### Run dev infrastructure with docker.
You can do it with or without the ELK stack.

To run is without (default):
```
cd ./infrastructure
docker-compose build
sh up-wo-elk.sh
```

To run with ELK, enable `logstash` in nomad.config:logstash, and start the docker compose with
```
docker-compose up
```
You can reach the Kibana with [localhost:5601](http://localhost:5601).
The index prefix for logs is `logstash-`.

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
