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

# This dockerfile describes an image that can be used to run the
# - nomad processing worker
# - nomad upload handler that initiates processing after upload
# - nomad api

# The dockerfile is multistaged to use a fat, more convinient build image and
# copy only necessities to a slim final image

# We use slim for the final image
FROM python:3.6-slim as final

# First, build everything in a build image
FROM python:3.6-stretch as build
# Make will be necessary to build the docs with sphynx
RUN apt-get update && apt-get install -y make
RUN mkdir /install
WORKDIR /install

# We also install the -dev dependencies, to use this image for test and qa
RUN pip install --upgrade pip
COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt

# Use docker build --build-args CACHEBUST=2 to not cache this (e.g. when you know deps have changed)
ARG CACHEBUST=1

# Install all NOMAD-CoE dependencies and nomad@FAIRDI
COPY . /install
RUN sh dependencies.sh
RUN pip install .
WORKDIR /install/docs
RUN make html
RUN \
    find /usr/local/lib/python3.6/ -name 'tests' ! -path '*/networkx/*' -exec rm -r '{}' + && \
    find /usr/local/lib/python3.6/ -name 'test' -exec rm -r '{}' + && \
    find /usr/local/lib/python3.6/site-packages/ -name '*.so' -print -exec sh -c 'file "{}" | grep -q "not stripped" && strip -s "{}"' \;

# Second, create a slim final image
FROM final

RUN apt-get update && apt-get install -y --no-install-recommends libgomp1 && apt-get install -y libmagic-dev

# copy the sources for tests, coverage, qa, etc.
COPY . /app
WORKDIR /app
# transfer installed packages from dependency stage
COPY --from=build /usr/local/lib/python3.6/site-packages /usr/local/lib/python3.6/site-packages
# copy the meta-info, since it files are loaded via relative paths. TODO that should change.
COPY --from=build /install/dependencies/nomad-meta-info /app/dependencies/nomad-meta-info
# copy the documentation, its files will be served by the API
COPY --from=build /install/docs/.build /app/docs/.build
# copy the nomad command
COPY --from=build /usr/local/bin/nomad /usr/local/bin/nomad

RUN mkdir -p /app/.volumes/fs
RUN useradd -ms /bin/bash nomad
RUN chown -R nomad /app
USER nomad

VOLUME /app/.volumes/fs

EXPOSE 8000
