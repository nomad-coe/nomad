import re
import unidecode

from nomad import infrastructure

infrastructure.setup_logging()

existing = set()

for user in infrastructure.keycloak.search_user(max=2000):
    if not re.match(r'^[a-zA-Z0-9_\-\.]+$', user.username):
        # need to replace username
        if user.first_name is not None and user.last_name is not None:
            user.username = '%s%s' % (user.first_name[:1], user.last_name)
        elif user.last_name is not None:
            user.username = user.last_name
        elif '@' in user.username:
            user.username = user.username.split('@')[0]

        user.username = unidecode.unidecode(user.username.lower())
        user.username = re.sub('[^0-9a-zA-Z_\-\.]+', '', user.username)

    index = 1
    while user.username in existing:
        user.username += '%d' % index
        index += 1

    existing.add(user.username)
    infrastructure.keycloak._admin_client.update_user(
        user_id=user.user_id, payload=dict(username=user.username))

    print(user.username)
