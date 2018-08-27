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
#
# The dockerfile is multistages to getaway with building on a larger base image.

# using the base image with most requirements already installed
FROM nomad-xt_requirements:latest as requirements
# we use slim for the final image
FROM python:3.6-slim as final


# dependency stage is used to install nomad coe projects
FROM requirements as dependencies

# do stuff
WORKDIR /install
COPY nomad/dependencies.py nomad/dependencies.py
COPY nomad/config.py nomad/config.py
RUN python nomad/dependencies.py


# last stage is used to install the actual code, nomad user, volumes
FROM final
# transfer installed packages from dependency stage
COPY --from=dependencies /usr/local/lib/python3.6/site-packages /usr/local/lib/python3.6/site-packages
# we also need to copy the install dir, since nomad coe deps are installed with -e
# TODO that should be changed in production!
COPY --from=dependencies /install /install

# do stuff
COPY . /app
WORKDIR /app
RUN pip install -e .
WORKDIR /app/docs
RUN make html
WORKDIR /app
RUN useradd -ms /bin/bash nomad
RUN mkdir -p /app/.volumes/fs; chown -R nomad /app/.volumes/fs
USER nomad
VOLUME /app/.volumes/fs
