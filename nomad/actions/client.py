import threading

import temporalio.converter
from temporalio.client import Client, TLSConfig
from temporalio.contrib.opentelemetry import OpenTelemetryPlugin
from temporalio.contrib.pydantic import PydanticPayloadConverter
from temporalio.runtime import PrometheusConfig, Runtime, TelemetryConfig

from nomad.actions._codec import EncryptionCodec
from nomad.config import config
from nomad.config.models.config import ModeEnum

_runtime: Runtime | None = None
_runtime_lock = threading.Lock()


def _get_metrics_config() -> PrometheusConfig | None:
    bind_address = config.temporal.prometheus_bind_address
    if bind_address is None:
        return None

    # The Temporal Prometheus exporter binds a process-local HTTP server.
    # The web app may run multiple Gunicorn workers in one deployment, so the
    # dedicated Temporal workers are the safer place to expose this endpoint.
    if config.meta.service == 'app':
        return None

    return PrometheusConfig(bind_address=bind_address)


def _get_runtime() -> Runtime:
    global _runtime

    if _runtime is None:
        with _runtime_lock:
            if _runtime is None:
                _runtime = Runtime(
                    telemetry=TelemetryConfig(metrics=_get_metrics_config()),
                    worker_heartbeat_interval=None,
                )

    return _runtime


def _load_cert_or_key(val: str | None) -> bytes | None:
    if not val:
        return None
    import os

    if os.path.exists(val):
        with open(val, 'rb') as f:
            return f.read()
    return val.encode('utf-8')


async def get_client() -> Client:
    # Ensure telemetry is initialized in this process before Temporal plugin checks
    # the global tracer provider type.
    if config.telemetry.enabled:
        from nomad.tracing import setup_tracing

        setup_tracing()

    host = f'{config.temporal.host}:{config.temporal.port}'
    data_converter = temporalio.converter.DataConverter(
        payload_converter_class=PydanticPayloadConverter,
        payload_codec=None
        # Disable encryption in dev mode
        if config.services.mode == ModeEnum.DEVELOPMENT
        else EncryptionCodec(),
    )
    plugins = []
    if config.telemetry.enabled:
        plugins.append(OpenTelemetryPlugin(add_temporal_spans=True))

    tls: bool | TLSConfig | None = None
    if (
        config.temporal.use_tls
        or config.temporal.tls_client_cert
        or config.temporal.tls_client_key
        or config.temporal.tls_server_root_ca_cert
        or config.temporal.tls_domain
    ):
        if (
            config.temporal.tls_client_cert
            or config.temporal.tls_client_key
            or config.temporal.tls_server_root_ca_cert
            or config.temporal.tls_domain
        ):
            tls = TLSConfig(
                client_cert=_load_cert_or_key(config.temporal.tls_client_cert),
                client_private_key=_load_cert_or_key(config.temporal.tls_client_key),
                server_root_ca_cert=_load_cert_or_key(
                    config.temporal.tls_server_root_ca_cert
                ),
                domain=config.temporal.tls_domain,
            )
        else:
            tls = True

    client = await Client.connect(
        host,
        namespace=config.temporal.namespace,
        api_key=config.temporal.api_key,
        tls=tls,
        data_converter=data_converter,
        runtime=_get_runtime(),
        plugins=plugins,
    )
    return client
