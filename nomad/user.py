import sys
import time
from mongoengine import Document, EmailField, StringField, ReferenceField, ListField
from passlib.apps import custom_app_context as pwd_context
from itsdangerous import TimedJSONWebSignatureSerializer as Serializer, BadSignature, SignatureExpired

from nomad import config


class User(Document):
    """ Represents users in the database. """
    email = EmailField(primary_key=True)
    name = StringField()
    password_hash = StringField()

    def hash_password(self, password):
        self.password_hash = pwd_context.encrypt(password)

    def verify_password(self, password):
        return pwd_context.verify(password, self.password_hash)

    def generate_auth_token(self, expiration=600):
        s = Serializer(config.services.api_secret, expires_in=expiration)
        return s.dumps({'id': self.id})

    @staticmethod
    def verify_auth_token(token):
        s = Serializer(config.services.api_secret)
        try:
            data = s.loads(token)
        except SignatureExpired:
            return None    # valid token, but expired
        except BadSignature:
            return None    # invalid token

        return User.objects(email=data['id']).first()


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

# provid a test user for testing
me = None


def ensure_test_users():
    global me
    me = User.objects(email='me@gmail.com').first()
    if me is None:
        me = User(
            email='me@gmail.com',
            name='Me Meyer')
        me.hash_password('nomad')
        me.save()
        time.sleep(1)

if 'sphinx' not in sys.modules:
    ensure_test_users()
