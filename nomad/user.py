import sys
import time
from mongoengine import Document, EmailField, StringField, ReferenceField, ListField


class User(Document):
    """ Represents users in the database. """
    email = EmailField(primary_key=True)
    name = StringField()


class DataSet(Document):
    name = StringField()
    description = StringField()
    doi = StringField()

    user = ReferenceField(User)
    calcs = ListField(StringField)

    meta = {
        'indexes': [
            'user',
            'doi',
            'calcs'
        ]
    }

# provid a fake user for testing
me = None


def ensure_test_users():
    global me
    me = User.objects(email='me@gmail.com').first()
    if me is None:
        me = User(email='me@gmail.com', name='Me Meyer')
        me.save()
        time.sleep(1)

if 'sphinx' not in sys.modules:
    ensure_test_users()
