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

### Build and run the dev infrastructure with docker.
There are different modes to work and develop with nomad-xt. First, you
can run everything in docker containers. Second, you can only run the 3rd party
services, like databases and message queues in docker, and run everything else (partially)
manually.

First, you need to build an image with all dependencies. We separated this so we
do not need to redownload and bwheel-build rarely changing dependencies all the time (and
multi stage builds and docker caching seems not to be good enough here).
From the root:
```
docker build -t nomad-requirements -f requirements.Dockerfile .
```

Now we can build the docker compose that contains all external services (redis, rabbitmq,
mongo, elastic, minio, elk) and nomad services (worker, handler, api, gui).
```
cd ./infrastructure
docker-compose build
```

You can run all containers with:
```
docker-compose up
```

You can alos run services selectively, e.g.
```
docker-compose up -d redis, rabbitmq, minio, minio-config, mongo, elastic, elk
docker-compose up worker handler
docker-compose up api gui
```

If you run the ELK stack (and enable logstash in nomad/config.py),
you can reach the Kibana with [localhost:5601](http://localhost:5601).
The index prefix for logs is `logstash-`.

If you want to access the minio object storage via the mc client, register the
infrastructure's minio host to the minio client (mc).
```
mc config host add minio http://localhost:9007 AKIAIOSFODNN7EXAMPLE wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY
```

### Serving minio, api, gui from one nginx
The docker-compose is setup to server all client accessible services from one webserver
via nginx *proxy_pass* directives.
This is currelty part of the gui image/container. The api is served at `/nomadxt/api`,
minio ad `/nomadxt/objects` and the gui ad `/nomadxt`. This also means that there
is some URL rewriting and prefixing in the api and gui.

### Run the nomad worker manually
You can run the worker, handler, api, and gui as part of the docker infrastructure, like
seen above.

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

### Run the handler
The handler is a small deamon that takes minio events and initiates the processing.
(This should actually be replaces, by minio talking to rabbitmq directly.)

You cna run the *handler* from docker-compose
```
docker-compose up nomad-handler
```

Or manually form the root
```
python -m nomad.handler
```

### Run the api
Either with docker, or:
```
python nomad/api.py
```

### Run the gui
When you run the gui on its own (e.g. with react dev server below), you have to have
the API running manually also. This *inside docker* API is configured for nging paths
and proxies, which are run by the gui container. But you can run the *production* gui
in docker and the dev server gui in parallel with an API in docker.
Either with docker, or:
```
cd gui
yarn
yarn start
```

### Run tests.
You need to have the infrastructure partially running: minio, elastic, rabbitmq, redis.
The rest should be mocked or provied by the tests. Make sure that you do no run any
worker and handler in parallel, as they will fight for tasks in the queue.
```
cd instrastructure
docker-compose up -d minio elastic rabbitmq, redis
cd ..
pytest -sv tests
```
