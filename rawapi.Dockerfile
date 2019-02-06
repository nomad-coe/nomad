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

# This dockerfile creates a limited api container that only runs the raw file api endpoint

# We use slim for the final image
FROM python:3.6-slim as final

# First, build everything in a build image
FROM python:3.6-stretch as build
RUN mkdir /install
WORKDIR /install

# We also install the -dev dependencies, to use this image for test and qa
RUN pip install --upgrade pip
COPY rawapi.requirements.txt rawapi.requirements.txt
RUN pip install -r rawapi.requirements.txt

# do that after the dependencies to use docker's layer caching
COPY . /install
RUN echo "from .app import app\nfrom . import raw" > /install/nomad/api/__init__.py
RUN pip install .
RUN \
    find /usr/local/lib/python3.6/ -name 'tests' ! -path '*/networkx/*' -exec rm -r '{}' + && \
    find /usr/local/lib/python3.6/ -name 'test' -exec rm -r '{}' + && \
    find /usr/local/lib/python3.6/site-packages/ -name '*.so' -print -exec sh -c 'file "{}" | grep -q "not stripped" && strip -s "{}"' \;

# Second, create a slim final image
FROM final
# copy the sources for tests, coverage, qa, etc.
COPY --from=build /install /app
WORKDIR /app
# transfer installed packages from dependency stage
COPY --from=build /usr/local/lib/python3.6/site-packages /usr/local/lib/python3.6/site-packages

RUN mkdir -p /raw
RUN useradd -ms /bin/bash nomad
RUN chown -R nomad /app
USER nomad

ENV NOMAD_FILES_OBJECTS_DIR /raw
ENV NOMAD_FILES_RAW_BUCKET data
ENV NOMAD_SERVICE rawapi

CMD python -m gunicorn.app.wsgiapp -b 0.0.0.0:8000 nomad.api:app

VOLUME /raw

EXPOSE 8000
