import pytest
from aiohttp import ClientSession

from nomad.actions.workers.health_check import (
    should_start_health_server,
    start_health_server,
)
from nomad.config import config
from nomad.config.models.config import ModeEnum, WorkerConfig


@pytest.mark.asyncio
async def test_start_health_server_serves_health(unused_tcp_port):
    runner = await start_health_server('127.0.0.1', unused_tcp_port)

    try:
        async with ClientSession() as session:
            async with session.get(
                f'http://127.0.0.1:{unused_tcp_port}/health'
            ) as response:
                assert response.status == 200
                assert await response.text() == 'OK'
    finally:
        await runner.cleanup()


def test_should_start_health_server_respects_development_mode(monkeypatch):
    worker_config = WorkerConfig(healthcheck_enabled=True)

    monkeypatch.setattr(config.services, 'mode', ModeEnum.PRODUCTION)
    assert should_start_health_server(worker_config)

    monkeypatch.setattr(config.services, 'mode', ModeEnum.DEVELOPMENT)
    assert not should_start_health_server(worker_config)


def test_should_start_health_server_respects_worker_config(monkeypatch):
    monkeypatch.setattr(config.services, 'mode', ModeEnum.PRODUCTION)

    assert not should_start_health_server(WorkerConfig(healthcheck_enabled=False))
