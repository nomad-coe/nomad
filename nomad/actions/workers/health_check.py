from aiohttp import web

from nomad.config import config
from nomad.config.models.config import ModeEnum, WorkerConfig


def should_start_health_server(worker_config: WorkerConfig) -> bool:
    return (
        worker_config.healthcheck_enabled
        and config.services.mode != ModeEnum.DEVELOPMENT
    )


# There's no built in health check in Temporal, but we can use the following endpoint to setup
# a healthcheck and then Docker/Kubernetes would detect worker failures/crashes and restart the container/pod.
# This approach is suggested in https://temporal.io/blog/deploying-temporal-workers-to-amazon-ecs
async def start_health_server(host: str, port: int) -> web.AppRunner:
    async def handle_health(_):
        return web.Response(text='OK')

    app = web.Application()
    app.router.add_get('/health', handle_health)

    runner = web.AppRunner(app)
    await runner.setup()

    site = web.TCPSite(runner, host, port)
    await site.start()

    return runner
