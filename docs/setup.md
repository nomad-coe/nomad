# Setup

### Install intra nomad dependencies.
This includes parsers, normalizers, python-common, meta-info, etc.
Those dependencies are managed and configures via python scripts.

This step is some what optional. Those dependencies are only needed for processing.
If you do not develop on the processing and do not need to run the workers from
your environment, and only use the docker image for processing, you can skip.

Install some pre-requisite requriements
```
pip install -r requirements-dep.txt
```

Run the dependency installation
```
python nomad/dependencies.py
```

### Install the the actual code

```
pip install -r requirements.txt
pip install -e .
```

### Run dev infrastructure with docker.
First, build the docker images
```
cd ./infrastructure
docker-compose build
```

This will download images for services like redis, minio, rabbitmq. It will configure
an existing image for the ELK stack. It will build the processing image that contains
all intra nomad dependencies (parsers, normalizers, etc.) and will run the celery workers
that do the processing.

You can run all containers, including ELK and processing workers:
```
docker-compose up
```

You can alos run services selectively, e.g.
```
docker-compose up redis, rabbitmq, minio
```

If you run the ELK stack (and enable logstash in nomad/config.py),
you can reach the Kibana with [localhost:5601](http://localhost:5601).
The index prefix for logs is `logstash-`.

If you want to access the minio object storage via the mc client, register the
infrastructure's minio host to the minio client (mc).
```
mc config host add minio http://localhost:9007 AKIAIOSFODNN7EXAMPLE wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY
```

### Run the celery worker
You can run the worker as part of the docker infrastructure.
```
cd infrastructure
docker-compose up nomad-worker
```

You can also run the worker yourself, e.g. to develop on the processing. To simply
run a worker do (from the root)
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

### Run tests.
You need to have the infrastructure running (including the nomad-worker service)
```
cd instrastructure
docker-compose up -d
cd ..
python tests/test_files.py
```
