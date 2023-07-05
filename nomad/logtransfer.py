import logging
import socketserver
from logging.handlers import RotatingFileHandler
import os
import re
import threading
import requests
import io
import gzip
from datetime import datetime, timedelta

from nomad import config

# It is important to make sure that logstash itself is disabled (otherwise it can happen that logs from the logtransfer
# server are also included to the logs
config.logstash.enabled = False


# logger only for the logstash server (none of the logs should be transferred to the central Nomad instance)
server_logger = logging.getLogger(name='logtransfer_server')
server_logger.setLevel(logging.INFO)
server_handler = logging.StreamHandler()
server_handler.setFormatter(logging.Formatter('%(levelname)s: %(asctime)s %(message)s', '%Y-%m-%d %H:%M:%S'))
server_logger.addHandler(server_handler)


def get_log_filepath():
    configured_path = os.path.join(config.fs.tmp, config.logtransfer.log_filename)
    return os.path.abspath(configured_path)


def get_all_rotated_logfiles():
    '''
    Returns all rotated logfiles.

    A rotated logfile has the following structure:

    [logfilename].txt.[index]

    where [index] is an integer (the higher the index the older the logs in the file).

    In contrast, for the active logfile, where logs may still be written, [index] is empty.
    The list is sorted such that the oldest file with the highest index appears first.
    '''
    rotated_logfiles = []

    for filepath in os.listdir(config.fs.tmp):
        if re.search(fr'{config.logtransfer.log_filename}.[0-9]+', filepath):
            fullpath = os.path.abspath(os.path.join(config.fs.tmp, filepath))
            n_rotated_logfiles = int(fullpath.split('.')[-1])
            rotated_logfiles.append((fullpath, n_rotated_logfiles))

    return [f[0] for f in sorted(rotated_logfiles, key=lambda x: x[1], reverse=True)]


def gzip_bytes(msg) -> bytes:
    buf = io.BytesIO()
    with gzip.GzipFile(mode='wb', fileobj=buf) as f:
        f.write(msg)
    return buf.getvalue()


def clear_logfiles():
    '''
    Perform a rollover on the current logfile and then remove all rotated logfiles.

    Note that if the logger is still used than this operation is unsafe and potentially breaks the logger. This
    function is mainly intended for testing.
    '''

    try:
        logfile = get_log_filepath()
        os.remove(logfile)
    except FileNotFoundError:
        pass  # don't worry

    for f in get_all_rotated_logfiles():
        try:
            os.remove(f)
        except FileNotFoundError:
            pass  # don't worry


def is_empty_logfile():
    logfile = get_log_filepath()
    if not os.path.exists(logfile):
        return True
    else:
        return os.stat(logfile).st_size == 0


class NotifyRotatingFileHandler(RotatingFileHandler):
    '''Adapted `RotatingFileHandler` used within in the `Logtransfer` class.

    This adaptation of RotatingFileHandler does not perform a file rollover when `maxBytes` is reached. Instead, it
    notifies the Logtransfer class thread. This is because data is being sent to the central Nomad, and performing a
    file rollover at the same time could result in file conflicts and race conditions.

    The Lock is needed to protect `shouldRollover()` and `doRollover()` as they are not thread-safe (see attribute
    `stream` in the handler).'''
    rollover_lock = threading.Lock()
    notify_rollover_event: threading.Event = None

    def set_event_attribute(self, event: threading.Event):
        '''Set the event to the class on which this handler can notify the thread to safely perform the file
        rollover.'''
        self.notify_rollover_event = event

    def is_event_attribute_set(self):
        return isinstance(self.notify_rollover_event, threading.Event)

    def emit(self, record) -> None:
        '''Emit a record.

        Output the record to the file, if `maxBytes` (see `RotatingFileHandler`) is reached a
        `threading.Event` is set (for the event see `Logtransfer.transfer_logs`.
        '''
        try:
            with self.rollover_lock:  # shouldRollover and doRollover are not thread-safe
                # mypy does not recognize the method in the superclass
                if self.shouldRollover(record):  # type: ignore
                    # print(f'notify main logstransfer thread self.stream.tell()={self.stream.tell()}')
                    # notify a thread which is responsible for the rollover (and other tasks)
                    self.notify_rollover_event.set()

                # write to file
                logging.FileHandler.emit(self, record)
        except Exception:
            self.handleError(record)


class Logtransfer:
    '''
    Responsible for rotating the logfile and transfering the logs to central Nomad.

    The main method is `transfer_logs` which is intended to run in a thread (concurrently to the logstash proxy server
    that receives and writes new logs).
    '''

    def __init__(
            self,
            rotating_file: NotifyRotatingFileHandler):

        self.reached_max_bytes_event = threading.Event()
        self.rotating_file = rotating_file
        self.rotating_file.set_event_attribute(self.reached_max_bytes_event)

        self.submit_interval = config.logtransfer.submit_interval

        if self.submit_interval <= 0 and not isinstance(self.submit_interval, int):
            raise TypeError(
                f'submit_interval must be a positive integer value. Got {self.submit_interval} '
                f'with type {type(self.submit_interval)}')

        self.central_nomad_api_url = config.oasis.central_nomad_deployment_url

        if not self.central_nomad_api_url.endswith('/'):
            self.central_nomad_api_url += '/'

        self.raise_unexpected_exceptions = config.logtransfer.raise_unexpected_exceptions

    def _rotated_logfile_iterator(self):

        all_rotated_logfiles = get_all_rotated_logfiles()

        if len(all_rotated_logfiles) > 0:
            server_logger.info(
                f'collected files: {[os.path.basename(f) for f in all_rotated_logfiles]} in '
                f'directory {config.fs.tmp}'
            )

            for path in all_rotated_logfiles:
                with open(path, 'rb') as f:
                    file_content = f.read()

                yield file_content, path
        else:
            server_logger.info('No logfiles to submit.')

    def _submit_logfile_to_central(self, logfile_content: bytes):

        central_post_federation_url = f'{self.central_nomad_api_url}federation/logs/'

        try:
            headers = {'Content-Encoding': 'gzip'}

            ret_request = requests.post(
                central_post_federation_url,
                headers=headers,
                data=gzip_bytes(logfile_content)
            )

            is_successful = ret_request.status_code == 200

            if is_successful:
                submitted_bytes = ret_request.json()['filesize']
                server_logger.info(
                    f'Successfully submitted logfile ({submitted_bytes} bytes) with HTTP status code '
                    f'{ret_request.status_code} to central Oasis at {ret_request.url}.'
                )
            else:
                server_logger.info(
                    f'Submission of logfiles to {central_post_federation_url} failed with HTTP '
                    f'status code {ret_request.status_code}. \n '
                    'logfiles will be included again in next submission.'
                )

        except requests.exceptions.ConnectionError:
            is_successful = False
            server_logger.info(f'HTTP connection to {central_post_federation_url} could not be established')

        return is_successful

    def _remove_logfile(self, filepath):
        try:
            os.remove(filepath)
        except FileNotFoundError:
            pass  # this should not happen but is also not really a problem

    def transfer_logs(self) -> None:
        # (mypy does not recognize the superclass attribute)
        server_logger.info(
            f'Start logtransfer thread with '  # type: ignore
            f'{self.submit_interval=} sec | '
            f'{self.rotating_file.maxBytes=} bytes'
        )

        while True:
            try:
                for file_content, filepath in self._rotated_logfile_iterator():
                    if file_content.strip() != b'':
                        is_successful = self._submit_logfile_to_central(file_content)

                        if is_successful:
                            self._remove_logfile(filepath)
                        else:
                            break  # do not attempt to send other files, try next time...
                    else:
                        # file seems to be empty
                        self._remove_logfile(filepath)
            except Exception as e:
                if not self.raise_unexpected_exceptions:
                    # do not kill the thread and hope it will fix for the next iteration
                    server_logger.info(f'unexpected exception was raised \n{e}')
                else:
                    raise e

            next_submission = datetime.now() + timedelta(seconds=self.submit_interval)
            next_submission = next_submission.strftime('%Y-%m-%dT%H:%M:%S')  # type: ignore
            server_logger.info(f'The next planned submission is at {next_submission}.')

            self.reached_max_bytes_event.clear()
            self.reached_max_bytes_event.wait(self.submit_interval)
            is_notified = self.reached_max_bytes_event.is_set()

            server_logger.info(f'Thread is waking up (is_notified={is_notified}).')

            with self.rotating_file.rollover_lock:
                # mypy does not recognize the superclass attribute
                if self.rotating_file.stream.tell() != 0:  # type: ignore
                    server_logger.info(
                        'Perform rollover on logfile. stream position in '  # type: ignore
                        f'bytes={self.rotating_file.stream.tell()}')

                    # Note: this should be the only place where doRollover() is called. A call outside this thread
                    # causes a race condition on the methods where the logfiles are submitted and removed
                    self.rotating_file.doRollover()
                else:
                    server_logger.info(
                        'No rollover on logfile performed because rotating_file.stream.tell() is at position 0.'
                    )


class LogstashTCPHandler(socketserver.StreamRequestHandler):

    def handle(self) -> None:
        self.server: TCPServerReuseAddress

        while True:
            logline = self.rfile.readline()

            if logline == b'':
                # "empty byte" signals that client closed connection, nothing more to read
                break

            # remove newlines from the right, because a newline is also internally included in the logstash_logger
            logline = logline.rstrip()

            if len(logline) > 0:
                # print(f'logtransfer received {logline}')
                self.server.logger.info(logline.decode())

        self.server.rotating_file_handler.flush()


class TCPServerReuseAddress(socketserver.TCPServer):
    '''
    This class overwrites default class parameters of TCPServer.
    This is recommended at https://stackoverflow.com/a/42147927)
    '''

    # TODO: consider also if it is better to inherit from ThreadingTCPServer because (from Python docu):
    #  https://docs.python.org/3/library/socketserver.html#socketserver.ThreadingTCPServer
    #  On the other hand, if you are building an HTTP server where all data is stored externally (for instance, in the
    #  file system), a synchronous class will essentially render the service “deaf” while one request is being handled
    #  – which may be for a very long time if a client is slow to receive all the data it has requested. Here a
    #  threading or forking server prevents this.
    #  -- BUT NOTE: this may not work for Windows, because the ThreadingMixIn is only supported for POSIX systems

    # this allows to quickly set up the server again after it may have terminated
    allow_reuse_port = True
    allow_reuse_address = True
    request_queue_size = 20  # default is 5, set to higher value here to avoid ConnectionRefused errors

    def __init__(self, logger, *args, **kwargs):
        super(TCPServerReuseAddress, self).__init__(*args, **kwargs)
        self.logger = logger
        self.rotating_file_handler = self.logger.handlers[0]


def _start_logstash_proxy_server(logger, host, port):

    if isinstance(port, str):
        # accept parseable str to align with Nomad config specification
        try:
            port = int(port)
        except ValueError:
            raise TypeError(
                f'Server port must be of type integer (or parseable string). Got port={port} of type {type(port)}'
            )

    logstash_proxy_server = TCPServerReuseAddress(logger, (host, port), LogstashTCPHandler, bind_and_activate=True)
    logstash_proxy_server.timeout = 30  # is closed after timeout period with no requests being received

    with logstash_proxy_server as server:
        server_logger.info(f'Start logstash proxy server on host={host} port={port}.')

        try:
            # NOTE: don't use the serve_forever() method as it does not account for the timeout
            while True:
                server.handle_request()
        except Exception as e:
            server_logger.info('logstash proxy server shutdown unexpectedly')
            raise e


def _initialize_logger_and_handler():
    server_logger.info(
        f'Initialize logger and handler. Location of intermediate logfiles before '
        f'submission:\n{get_log_filepath()}'
    )
    logstash_logger = logging.getLogger(name='logstash')
    # Note that backupCount must be a positive number (otherwise no file rollovers is performed).
    # the number should for performance reasons not be too large. If backupCount is reached then the oldest log file
    # is deleted. See https://stackoverflow.com/a/56954614
    rotating_logfile_handler = NotifyRotatingFileHandler(
        filename=get_log_filepath(),
        mode='a',
        backupCount=config.logtransfer.backup_count,
        maxBytes=config.logtransfer.max_bytes
    )

    logstash_logger.addHandler(rotating_logfile_handler)
    rotating_logfile_handler.setLevel(logging.INFO)
    logstash_logger.setLevel(logging.INFO)

    return logstash_logger, rotating_logfile_handler


def start_logtransfer_service(host=None, port=None):
    '''
    Start logtransfer service.

    The two parameters host and port are mainly necessary for testing. If left None, then the values from the
    Nomad config are used.
    '''
    # for appropriate arguments see attributes in
    # nomad.config.logstash and nomad.config.logstash_proxy and config.oasis.central_nomad_api_url

    server_logger.info(f'Start logtransfer service')

    logstash_logger, rotating_logfile_handler = _initialize_logger_and_handler()

    # thread to frequently submit collected log data and submit them to central Nomad
    transfer_to_central = Logtransfer(rotating_file=rotating_logfile_handler)

    # is set within Logtransfer
    assert rotating_logfile_handler.is_event_attribute_set()

    d = threading.Thread(target=transfer_to_central.transfer_logs, name='logtransfer')
    d.setDaemon(True)  # kill thread immediately if main thread (running the server) terminates
    d.start()

    if host is None:
        host = config.logstash.host

    if port is None:
        port = config.logstash.tcp_port

    # start logstash proxy server to receive logs on main thread
    _start_logstash_proxy_server(logger=logstash_logger, host=host, port=port)
