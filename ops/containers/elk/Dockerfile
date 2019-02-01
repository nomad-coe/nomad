FROM sebp/elk

RUN apt-get update -qq \
&& apt-get install -qqy elasticsearch-curator

ENV LOGSTASH_PATH_CONF /etc/logstash

ADD ./rotate_indices.curator.yml /etc/curator/rotate_indices.yml
ADD ./curator.cron.sh /etc/cron.daily/curator
RUN chmod a+x /etc/cron.daily/curator

ADD ./02-beats-input.conf ${LOGSTASH_PATH_CONF}/conf.d/02-beats-input.conf
ADD ./30-output.conf ${LOGSTASH_PATH_CONF}/conf.d/30-output.conf
ADD ./kibana.yml ${KIBANA_HOME}/config/kibana.yml