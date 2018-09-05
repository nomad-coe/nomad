import sys
from mongoengine import Document, EmailField, StringField, ReferenceField, ListField


class User(Document):
    """ Represents users in the database. """
    email = EmailField(primary=True)
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
if 'sphinx' not in sys.modules:
    me = User.objects(email='me@gmail.com').first()
    if me is None:
        me = User(email='me@gmail.com', name='Me Meyer')
        me.save()
