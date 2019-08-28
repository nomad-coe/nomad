FROM jboss/keycloak
ARG admin_password=password
ARG prefix='fairdi\/keycloak'

USER jboss

RUN echo "s/<web-context>auth<\\/web-context>/<web-context>$prefix\\/auth<\\/web-context>/"
RUN sed -i -e "s/<web-context>auth<\\/web-context>/<web-context>$prefix\\/auth<\\/web-context>/" $JBOSS_HOME/standalone/configuration/standalone.xml
RUN sed -i -e "s/<web-context>auth<\\/web-context>/<web-context>$prefix\\/auth<\\/web-context>/" $JBOSS_HOME/standalone/configuration/standalone-ha.xml
RUN sed -i -e "s/name=\"\\/\"/name=\"\\/$prefix\\/\"/" $JBOSS_HOME/standalone/configuration/standalone.xml
RUN sed -i -e "s/name=\"\\/\"/name=\"\\/$prefix\\/\"/" $JBOSS_HOME/standalone/configuration/standalone-ha.xml
RUN sed -i -e "s/\\/auth/\\/$prefix\\/auth/" $JBOSS_HOME/welcome-content/index.html

ENV KEYCLOAK_USER=admin
ENV KEYCLOAK_PASSWORD=$admin_password
ENV KEYCLOAK_IMPORT=/opt/jboss/keycloak/fairdi_nomad_test.json
ENV PROXY_ADDRESS_FORWARDING=true

ADD ./fairdi_nomad_test.json /opt/jboss/keycloak/fairdi_nomad_test.json
ADD ./fairdi_nomad_prod.json /opt/jboss/keycloak/fairdi_nomad_prod.json

ADD ./material_theme /opt/jboss/keycloak/themes/material

RUN curl http://central.maven.org/maven2/org/mindrot/jbcrypt/0.4/jbcrypt-0.4.jar > /tmp/jbcrypt-0.4.jar
RUN /opt/jboss/keycloak/bin/jboss-cli.sh --command="module add --name=org.mindrot.jbcrypt --resources=/tmp/jbcrypt-0.4.jar"
RUN curl -L https://github.com/leroyguillaume/keycloak-bcrypt/releases/download/1.1.0/keycloak-bcrypt-1.1.0.jar > /opt/jboss/keycloak/standalone/deployments/keycloak-bcrypt-1.1.0.jar