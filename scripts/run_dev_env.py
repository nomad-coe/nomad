import sys

from nomad.cli.admin.run import run_appworker

# run both app and worker in development mode
# that is, use one worker for fastapi and another for celery
# to minimize the resources used
# to use custom host and port, run the following
#
# python run_dev_env.py --host=1.2.3.4 --port=5678
#

if __name__ == '__main__':
    host = None
    port = None

    for arg in sys.argv:
        if '--host' in arg:
            host = arg.split('=')[1]
        if '--port' in arg:
            port = arg.split('=')[1]

    run_appworker(dev=True, app_host=host, app_port=port)
