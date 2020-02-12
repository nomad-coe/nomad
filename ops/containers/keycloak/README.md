This basically uses the "official" 7.0.0 keycloak with some additions

- material theme
- bcrypt support (https://github.com/leroyguillaume/keycloak-bcrypt/)
- a custom form action to check the username (that has to be build before building the image!)
- realm definitions for prod and test realm
- scripts for import and export of all users

To build the image:
```
cd registration_form_action
mvn package
cd ..
docker build -t nomad/keycloak .
```