FROM sebp/elk

ENV LOGSTASH_PATH_CONF /etc/logstash

ADD ./02-beats-input.conf ${LOGSTASH_PATH_CONF}/conf.d/02-beats-input.conf
ADD ./30-output.conf ${LOGSTASH_PATH_CONF}/conf.d/30-output.conf
ADD ./kibana.yml ${KIBANA_HOME}/config/kibana.yml