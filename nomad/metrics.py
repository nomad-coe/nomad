#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import os
import time

from fastapi import FastAPI, Response

from nomad.config import config

# Set multiprocess directory before importing prometheus_client
if config.telemetry.metrics.api_prometheus_enabled:
    multiproc_dir = os.path.join(config.fs.tmp, 'prometheus_multiproc')
    os.makedirs(multiproc_dir, exist_ok=True)
    os.environ['PROMETHEUS_MULTIPROC_DIR'] = multiproc_dir

from prometheus_client import (
    CONTENT_TYPE_LATEST,
    CollectorRegistry,
    Counter,
    Gauge,
    Histogram,
    generate_latest,
    multiprocess,
)

REQUEST_COUNT = Counter(
    'nomad_fastapi_requests_total',
    'Total count of FastAPI requests.',
    ['method', 'path', 'status_code'],
)

REQUEST_LATENCY = Histogram(
    'nomad_fastapi_request_duration_seconds',
    'Request latency in seconds.',
    ['method', 'path'],
    buckets=(
        0.005,
        0.01,
        0.025,
        0.05,
        0.075,
        0.1,
        0.25,
        0.5,
        0.75,
        1.0,
        2.5,
        5.0,
        7.5,
        10.0,
        float('inf'),
    ),
)

REQUEST_SIZE = Histogram(
    'nomad_fastapi_request_size_bytes',
    'Size of incoming FastAPI request bodies in bytes.',
    ['method', 'path'],
    buckets=(
        128,
        512,
        2048,
        8192,
        32768,
        131072,
        524288,
        2097152,
        8388608,
        float('inf'),
    ),
)

RESPONSE_SIZE = Histogram(
    'nomad_fastapi_response_size_bytes',
    'Size of outgoing FastAPI response bodies in bytes.',
    ['method', 'path'],
    buckets=(
        128,
        512,
        2048,
        8192,
        32768,
        131072,
        524288,
        2097152,
        8388608,
        float('inf'),
    ),
)

REQUESTS_IN_PROGRESS = Gauge(
    'nomad_fastapi_requests_in_progress',
    'Number of concurrent FastAPI requests in progress.',
    ['method'],
    multiprocess_mode='livesum',
)


class PrometheusASGIMiddleware:
    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        if scope['type'] != 'http':
            await self.app(scope, receive, send)
            return

        # Do not instrument the metrics endpoint itself
        path = scope.get('path', '')
        if path.endswith('/metrics'):
            await self.app(scope, receive, send)
            return

        # Extract request size from content-length header
        request_size = 0
        for name, value in scope.get('headers', []):
            if name == b'content-length':
                try:
                    request_size = int(value.decode('latin1'))
                except ValueError:
                    pass
                break

        REQUESTS_IN_PROGRESS.labels(method=scope['method']).inc()

        start_time = time.perf_counter()
        status_code = [500]
        response_size = [0]

        async def wrapped_send(message):
            if message['type'] == 'http.response.start':
                status_code[0] = message['status']
            elif message['type'] == 'http.response.body':
                body_chunk = message.get('body', b'')
                response_size[0] += len(body_chunk)
            await send(message)

        try:
            await self.app(scope, receive, wrapped_send)
        except Exception as exc:
            status_code[0] = 500
            raise exc
        finally:
            latency = time.perf_counter() - start_time

            # Resolve templated path to prevent cardinality explosion
            route = scope.get('route')
            if route and hasattr(route, 'path'):
                root_path = scope.get('root_path', '')
                templated_path = f'{root_path}{route.path}'
            else:
                templated_path = '/unmatched'

            REQUEST_COUNT.labels(
                method=scope['method'],
                path=templated_path,
                status_code=str(status_code[0]),
            ).inc()

            REQUEST_LATENCY.labels(method=scope['method'], path=templated_path).observe(
                latency
            )

            REQUEST_SIZE.labels(method=scope['method'], path=templated_path).observe(
                request_size
            )

            RESPONSE_SIZE.labels(method=scope['method'], path=templated_path).observe(
                response_size[0]
            )

            REQUESTS_IN_PROGRESS.labels(method=scope['method']).dec()


def setup_prometheus(app: FastAPI):
    if not config.telemetry.metrics.api_prometheus_enabled:
        return

    # Add raw ASGI middleware
    app.add_middleware(PrometheusASGIMiddleware)

    # Expose the /metrics endpoint
    app_base = config.services.api_base_path

    @app.get(f'{app_base}/metrics', include_in_schema=False)
    async def metrics():
        if 'PROMETHEUS_MULTIPROC_DIR' in os.environ:
            registry = CollectorRegistry()
            multiprocess.MultiProcessCollector(registry)
            data = generate_latest(registry)
        else:
            data = generate_latest()
        return Response(content=data, media_type=CONTENT_TYPE_LATEST)
