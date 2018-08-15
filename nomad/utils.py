from typing import Union, IO, cast
import hashlib
import base64


def hash(obj: Union[IO, str]) -> str:
    """ First 28 character of an URL safe base 64 encoded sha512 digest. """
    hash = hashlib.sha512()
    if getattr(obj, 'read', None) is not None:
        for data in iter(lambda: cast(IO, obj).read(65536), b''):
            hash.update(data)
    elif isinstance(obj, str):
        hash.update(obj.encode('utf-8'))

    return base64.b64encode(hash.digest(), altchars=b'-_')[0:28].decode('utf-8')
