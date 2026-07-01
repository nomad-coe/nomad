import pytest
from temporalio.api.common.v1 import Payload

from nomad.actions._codec import EncryptionCodec
from nomad.config import config
from nomad.config.models.config import ModeEnum


@pytest.fixture
def codec_config(monkeypatch):
    monkeypatch.setattr(config.services, 'mode', ModeEnum.DEVELOPMENT)
    monkeypatch.setattr(
        config.services, 'api_secret', 'local-secret-that-is-at-least-32-bytes'
    )
    monkeypatch.setattr(config.temporal, 'payload_codec_key', None)
    monkeypatch.setattr(config.temporal, 'payload_codec_key_id', 'default')


@pytest.mark.asyncio
async def test_configured_codec_key_is_recorded_and_used(codec_config, monkeypatch):
    monkeypatch.setattr(
        config.temporal,
        'payload_codec_key',
        'federation-secret-that-is-at-least-32-bytes',
    )
    monkeypatch.setattr(config.temporal, 'payload_codec_key_id', 'federation-v1')
    codec = EncryptionCodec()
    payload = Payload(metadata={'encoding': b'json/plain'}, data=b'input')

    encoded = await codec.encode([payload])

    assert encoded[0].metadata['encryption-key-id'] == b'federation-v1'
    assert await codec.decode(encoded) == [payload]


@pytest.mark.asyncio
async def test_configured_codec_key_can_decode_default_payload(
    codec_config, monkeypatch
):
    payload = Payload(metadata={'encoding': b'json/plain'}, data=b'input')
    default_payload = (await EncryptionCodec().encode([payload]))[0]
    assert 'encryption-key-id' not in default_payload.metadata

    monkeypatch.setattr(
        config.temporal,
        'payload_codec_key',
        'federation-secret-that-is-at-least-32-bytes',
    )
    monkeypatch.setattr(config.temporal, 'payload_codec_key_id', 'federation-v1')

    assert await EncryptionCodec().decode([default_payload]) == [payload]


@pytest.mark.asyncio
async def test_unknown_codec_key_id_is_rejected(codec_config, monkeypatch):
    monkeypatch.setattr(
        config.temporal,
        'payload_codec_key',
        'federation-secret-that-is-at-least-32-bytes',
    )
    monkeypatch.setattr(config.temporal, 'payload_codec_key_id', 'federation-v1')
    codec = EncryptionCodec()
    encoded = await codec.encode([Payload(data=b'input')])
    encoded[0].metadata['encryption-key-id'] = b'federation-v2'

    with pytest.raises(ValueError, match="key ID 'federation-v2'"):
        await codec.decode(encoded)
