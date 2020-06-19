FROM bitnami/matomo:3

RUN sed -i -e 's/action_name=/a=b\&action_name=/g' /opt/bitnami/matomo/piwik.js
RUN sed -i -e 's/&send_image=0/\&/g' /opt/bitnami/matomo/piwik.js