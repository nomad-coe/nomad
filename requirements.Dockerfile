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
#
# This dockerfile creates a python image with most rarely changing requirements
# preinstalled, to speed up the build of the actual nomad-xt image.

# some dependencies with complex wheels do not build with -alpine, -slim
FROM python:3.6-stretch as builder

# do stuff
RUN mkdir /install
WORKDIR /install
COPY requirements.txt requirements.txt
COPY requirements-dep.txt requirements-dep.txt
RUN pip install -r requirements.txt
RUN pip install -r requirements-dep.txt
