from nomad.utils import get_logger
from nomad.logtransfer import transfer_logs

get_logger(__name__).info('logger initialized')
transfer_logs()
