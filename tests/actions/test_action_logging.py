import logging
import os

import structlog

from nomad.actions.action_logging import WorkflowRoutingHandler
from nomad.config import config


def test_interceptor_and_routing_handler(tmp_path, monkeypatch):
    """
    Test the automatic routing via contextvars and standard logging handler.
    """
    monkeypatch.setattr(config.fs, 'actions', str(tmp_path))

    root_logger = logging.getLogger()

    # We clear out existing handlers for the test to avoid duplicate outputs
    # and we add our interceptor handler
    old_handlers = root_logger.handlers[:]
    root_logger.handlers = []

    handler = WorkflowRoutingHandler()
    root_logger.addHandler(handler)
    root_logger.setLevel(logging.INFO)

    try:
        # Simulate activity context variable
        structlog.contextvars.clear_contextvars()
        structlog.contextvars.bind_contextvars(workflow_id='test-routed-workflow')

        # Now emit a generic log from root logger
        root_logger.info('This should magically route!')

        # Test another workflow
        structlog.contextvars.bind_contextvars(workflow_id='another-workflow')
        root_logger.error('Something went wrong over here')

        expected_log_dir = os.path.join(str(tmp_path), 'logs')
        file1 = os.path.join(expected_log_dir, 'test-routed-workflow.log')
        file2 = os.path.join(expected_log_dir, 'another-workflow.log')

        assert os.path.exists(file1)
        assert os.path.exists(file2)

        with open(file1) as f:
            content1 = f.read()
            assert 'This should magically route!' in content1
            assert 'Something went wrong' not in content1

        with open(file2) as f:
            content2 = f.read()
            assert 'Something went wrong over here' in content2
            assert 'magically route' not in content2
    finally:
        # Restore old state
        root_logger.handlers = old_handlers
        structlog.contextvars.clear_contextvars()
