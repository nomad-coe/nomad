FROM jboss/keycloak:7.0.0
ARG admin_password=password
ARG prefix='fairdi\/keycloak'

USER jboss

# Change the 'web-context' to allow reverse proxy behind custom path
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
ENV DB_VENDOR=h2

# Add exported real/client configurations
ADD ./fairdi_nomad_test.json /opt/jboss/keycloak/fairdi_nomad_test.json
ADD ./fairdi_nomad_prod.json /opt/jboss/keycloak/fairdi_nomad_prod.json

# Add a custom theme to mimic the rest of nomad's GUI
ADD ./material_theme /opt/jboss/keycloak/themes/material

# Add custom registration validation
ADD ./registration_form_action/target/form-extension-1.0-SNAPSHOT.jar /opt/jboss/keycloak/standalone/deployments/form-extension-1.0-SNAPSHOT.jar

# Additing a bcrypt password hashing algorithm to reuse old repo users
ADD ./jbcrypt-0.4.jar /tmp/jbcrypt-0.4.jar
RUN /opt/jboss/keycloak/bin/jboss-cli.sh --command="module add --name=org.mindrot.jbcrypt --resources=/tmp/jbcrypt-0.4.jar"
ADD ./keycloak-bcrypt-1.1.0.jar /opt/jboss/keycloak/standalone/deployments/keycloak-bcrypt-1.1.0.jar
