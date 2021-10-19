# build stage
FROM ghcr.io/linuxserver/baseimage-ubuntu:focal as nodebuilder

RUN \
    echo "**** install build deps ****" && \
    apt-get update && \
    apt-get install -y \
        gnupg && \
    curl -s https://deb.nodesource.com/gpgkey/nodesource.gpg.key | apt-key add - && \
    echo 'deb https://deb.nodesource.com/node_14.x focal main' \
        > /etc/apt/sources.list.d/nodesource.list && \
    apt-get update && \
    apt-get install -y \
        g++ \
        gcc \
        libpam0g-dev \
        make \
        nodejs

COPY north-guacamole-client /gclient

RUN \
    echo "**** install node modules ****" && \
    cd /gclient && \
    npm install

# runtime stage
FROM ghcr.io/linuxserver/webtop:ubuntu-openbox

RUN rm -rf /gclient
COPY --from=nodebuilder /gclient /gclient

# ports and volumes
EXPOSE 8888
VOLUME /config
