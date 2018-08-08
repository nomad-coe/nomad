FROM python:3.6-stretch
RUN rm /bin/sh && ln -s /bin/bash /bin/sh

ADD ./ /nomad
WORKDIR /nomad

RUN pip install -r requirements.txt
RUN pip install -e .
RUN pip install -r requirements-worker.txt
RUN python nomad/parsers.py

RUN useradd -ms /bin/bash nomad
USER nomad
