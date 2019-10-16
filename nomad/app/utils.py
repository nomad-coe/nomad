# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from werkzeug.exceptions import HTTPException
from flask_restplus import fields
from datetime import datetime
import pytz
import inspect

from nomad import utils


def with_logger(func):
    """
    Decorator for endpoint implementations that provides a pre configured logger and
    automatically logs errors on all 500 responses.
    """
    signature = inspect.signature(func)
    has_logger = 'logger' in signature.parameters
    wrapper_signature = signature.replace(parameters=tuple(
        param for param in signature.parameters.values()
        if param.name != 'logger'
    ))

    def wrapper(*args, **kwargs):
        if has_logger:
            args = inspect.getcallargs(wrapper, *args, **kwargs)
            logger_args = {
                k: v for k, v in args.items()
                if k in ['upload_id', 'calc_id']}
            logger = utils.get_logger(__name__, **logger_args)
            args.update(logger=logger)
        try:
            return func(**args)
        except HTTPException as e:
            if getattr(e, 'code', None) == 500:
                logger.error('Internal server error', exc_info=e)
            raise e
        except Exception as e:
            logger.error('Internal server error', exc_info=e)
            raise e

    wrapper.__signature__ = wrapper_signature
    return wrapper


class RFC3339DateTime(fields.DateTime):

    def format(self, value):
        if isinstance(value, datetime):
            return super().format(value.replace(tzinfo=pytz.utc))
        else:
            return str(value)


rfc3339DateTime = RFC3339DateTime()
