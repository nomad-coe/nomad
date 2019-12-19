docker exec -it nomad_keycloak keycloak/bin/standalone.sh \
    -Djboss.socket.binding.port-offset=200 \
    -Dkeycloak.migration.action=import \
    -Dkeycloak.migration.provider=singleFile \
    -Dkeycloak.migration.realmName=fairdi_nomad_prod \
    -Dkeycloak.migration.strategy=OVERWRITE_EXISTING \
    -Dkeycloak.migration.file=/export/fairdi_nomad_prod_latest.json
