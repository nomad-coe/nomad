import temporalio.converter
from temporalio.client import Client
from temporalio.contrib.opentelemetry import OpenTelemetryPlugin
from temporalio.contrib.pydantic import PydanticPayloadConverter
from temporalio.runtime import PrometheusConfig, Runtime, TelemetryConfig

from nomad.actions._codec import EncryptionCodec
from nomad.config import config
from nomad.config.models.config import ModeEnum


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
    runtime = Runtime(
        telemetry=TelemetryConfig(
            metrics=PrometheusConfig(
                bind_address=config.temporal.prometheus_bind_address
            )
            if config.temporal.prometheus_bind_address is not None
            else None
        ),
        worker_heartbeat_interval=None,
    )
    plugins = []
    if config.telemetry.enabled:
        plugins.append(OpenTelemetryPlugin(add_temporal_spans=True))

    client = await Client.connect(
        host,
        namespace=config.temporal.namespace,
        data_converter=data_converter,
        runtime=runtime,
        plugins=plugins,
    )
    return client
