import base64
from collections.abc import Iterable

from cryptography.fernet import Fernet
from temporalio.api.common.v1 import Payload
from temporalio.converter import PayloadCodec

from nomad.auth.tokens import check_api_secret
from nomad.config import config

_ENCRYPTED_ENCODING = b'binary/encrypted'
_KEY_ID_METADATA = 'encryption-key-id'


def _create_fernet(secret: str) -> Fernet:
    # Fernet requires a URL-safe, base64-encoded 32-byte key. Pad shorter
    # secrets (as used in tests) and truncate longer ones.
    secret_bytes = secret.encode()
    padded_key = b'\0' * max(32 - len(secret_bytes), 0) + secret_bytes
    return Fernet(base64.urlsafe_b64encode(padded_key[:32]))


class EncryptionCodec(PayloadCodec):
    """A PayloadCodec that encrypts/decrypts all Payloads."""

    def __init__(self) -> None:
        super().__init__()

        # Payloads using the default key retain the existing metadata format so
        # previously recorded workflow histories remain decryptable during replay.
        check_api_secret()
        self.default_fernet = _create_fernet(config.services.api_secret)

        codec_key = config.temporal.payload_codec_key
        self.key_id = config.temporal.payload_codec_key_id if codec_key else None
        self.fernet = _create_fernet(codec_key) if codec_key else self.default_fernet

    async def encode(self, payloads: Iterable[Payload]) -> list[Payload]:
        """Encrypt all payloads during encoding."""
        encoded_payloads = []
        for payload in payloads:
            metadata = {'encoding': _ENCRYPTED_ENCODING}
            if self.key_id is not None:
                metadata[_KEY_ID_METADATA] = self.key_id.encode()
            encoded_payloads.append(
                Payload(
                    metadata=metadata,
                    data=self.encrypt(payload.SerializeToString()),
                )
            )
        return encoded_payloads

    async def decode(self, payloads: Iterable[Payload]) -> list[Payload]:
        """Decode all payloads decrypting those with expected encoding."""
        ret: list[Payload] = []
        for p in payloads:
            # Ignore ones without our expected encoding
            if p.metadata.get('encoding') != _ENCRYPTED_ENCODING:
                ret.append(p)
                continue

            payload_key_id = p.metadata.get(_KEY_ID_METADATA)
            if payload_key_id is None:
                fernet = self.default_fernet
            else:
                key_id = payload_key_id.decode()
                if key_id != self.key_id:
                    raise ValueError(
                        f'No Temporal payload codec key configured for key ID {key_id!r}'
                    )
                fernet = self.fernet

            ret.append(Payload.FromString(fernet.decrypt(p.data)))
        return ret

    def encrypt(self, data: bytes) -> bytes:
        """Return data encrypted."""
        return self.fernet.encrypt(data)

    def decrypt(self, data: bytes) -> bytes:
        """Return data decrypted."""
        return self.fernet.decrypt(data)
