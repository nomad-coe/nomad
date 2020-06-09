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
# - nomad app (incl serving the gui)

# The dockerfile is multistaged to use a fat, more convinient build image and
# copy only necessities to a slim final image

# We use stretch instead of slim to allow compoilation of some dependencies. We
# do not bother to have an additional copy to a final slim image, because the space
# savings are minimal
FROM python:3.7-stretch as final

# First built the GUI in a gui build image
FROM node:latest as gui_build
RUN mkdir -p /app
WORKDIR /app
ENV PATH /app/node_modules/.bin:$PATH
COPY gui/package.json /app/package.json
COPY gui/yarn.lock /app/yarn.lock
RUN yarn
COPY gui /app
RUN yarn run build

# Second, build all python stuff in a python build image
FROM final
RUN mkdir /app

# Install linux package dependencies
RUN apt-get update
RUN apt-get install -y --no-install-recommends libgomp1
RUN apt-get install -y libmagic-dev curl vim make cmake swig libnetcdf-dev

# Install some specific dependencies necessary for the build process
RUN pip install --upgrade pip
RUN pip install fastentrypoints
RUN pip install pyyaml
RUN pip install numpy

# Install some specific dependencies to make use of docker layer caching
RUN pip install cython>=0.19
RUN pip install pandas
RUN pip install h5py
RUN pip install hjson
RUN pip install scipy
RUN pip install scikit-learn==0.20.2
RUN pip install ase==3.19.0
RUN pip install Pint
RUN pip install matid
RUN pip install mdtraj
RUN pip install mdanalysis

# Install pymolfile (required by some parsers)
RUN git clone -b nomad-fair https://gitlab.mpcdf.mpg.de/nomad-lab/pymolfile.git
WORKDIR /pymolfile/
RUN python3 setup.py install
RUN rm -rf /pymolfile

# Copy files and install nomad@FAIRDI
WORKDIR /app
COPY . /app
RUN python setup.py compile
RUN pip install .[all]
RUN python setup.py sdist
RUN cp dist/nomad-lab-*.tar.gz dist/nomad-lab.tar.gz

WORKDIR /app/docs
RUN make html

# Remove unnessesary files
WORKDIR /app
RUN \
    find /usr/local/lib/python3.7/ -name 'tests' ! -path '*/networkx/*' -exec rm -r '{}' + && \
    find /usr/local/lib/python3.7/ -name 'test' -exec rm -r '{}' + && \
    find /usr/local/lib/python3.7/site-packages/ -name '*.so' -print -exec sh -c 'file "{}" | grep -q "not stripped" && strip -s "{}"' \;

# Setup directories, users, rights
RUN mkdir -p /app/.volumes/fs
RUN useradd -ms /bin/bash nomad
RUN chown -R nomad /app
RUN chmod a+rx run.sh
USER nomad

VOLUME /app/.volumes/fs

EXPOSE 8000
