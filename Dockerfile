# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# This dockerfile describes an image that can be used to run the nomad processing
# worker than are capable of running all necessary celery task involving
# upload processing, parsing, nomalization, packaging analysis files, etc.

# This dockerfile was build with the following advide for small and incremental builds
# https://blog.realkinetic.com/building-minimal-docker-containers-for-python-applications-37d0272c52f3

# some dependencies with complex wheels do not build with -alpine, -slim
FROM python:3.6-stretch as base

# builder stage is used to install requirments and nomad dependencies
# (parsers, normalizers, python-common, meta-info, etc.)
FROM base as builder

RUN mkdir /install
WORKDIR /install

COPY requirements.txt requirements.txt
COPY requirements-dep.txt requirements-dep.txt
COPY nomad/dependencies.py nomad/dependencies.py
COPY nomad/config.py nomad/config.py

RUN pip install -r requirements.txt
RUN pip install -r requirements-dep.txt
RUN python nomad/dependencies.py

# second stage is used to install the actual code and run the celery worker as nomad user
FROM base
# transfer installs from builder stage
COPY --from=builder /install /install
COPY --from=builder /usr/local/lib/python3.6/site-packages /usr/local/lib/python3.6/site-packages
COPY . /app
WORKDIR /app

RUN pip install -e .

RUN useradd -ms /bin/bash nomad
RUN mkdir -p /app/.volumes/fs; chown -R nomad /app/.volumes/fs
USER nomad
VOLUME /app/.volumes/fs
