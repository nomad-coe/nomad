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

FROM alpine:3.6
RUN apk --no-cache --update-cache --virtual=.build-dependencies add g++ gfortran file binutils musl-dev python3-dev openblas-dev libstdc++ openblas make
RUN ln -s /usr/include/locale.h /usr/include/xlocale.h
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN python -m ensurepip
RUN ln -s /usr/bin/pip3 /usr/bin/pip

RUN mkdir /install
WORKDIR /install
# We also install the -dev dependencies, to use this image for test and qa
COPY requirements-dev.txt requirements-dev.txt
RUN pip install -r requirements-dev.txt
COPY requirements-dep.txt requirements-dep.txt
RUN pip install -r requirements-dep.txt
COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt
# Use docker build --build-args CACHEBUST=2 to not cache this (e.g. when you know deps have changed)
ARG CACHEBUST=1
COPY nomad/dependencies.py /install/nomad/dependencies.py
COPY nomad/config.py /install/nomad/config.py
RUN python nomad/dependencies.py
RUN ls -la .dependencies/parsers/vasp/
RUN ls -la .dependencies/parsers/vasp/vaspparser/
# do that after the dependencies to use docker's layer caching
COPY . /app
RUN pip app .
WORKDIR /app/docs
RUN make html
RUN \
    find /usr/lib/python3.6/ -name 'tests' -exec rm -r '{}' + && \
    find /usr/lib/python3.6/ -name 'test' -exec rm -r '{}' + && \
    find /usr/lib/python3.6/site-packages/ -name '*.so' -print -exec sh -c 'file "{}" | grep -q "not stripped" && strip -s "{}"' \;

RUN rm /usr/include/xlocale.h && \
    apk del .build-dependencies

RUN mkdir -p /app/.volumes/fs
RUN mkdir -p /nomad
RUN useradd -ms /bin/bash nomad
RUN chown -R nomad /app
RUN chown -R nomad /nomad
USER nomad

VOLUME /app/.volumes/fs
VOLUME /nomad
