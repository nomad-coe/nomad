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

from structlog import BoundLogger
from flask_restplus import fields
from datetime import datetime
import pytz
from contextlib import contextmanager


logger: BoundLogger = None
''' A logger pre configured with information about the current request. '''


class RFC3339DateTime(fields.DateTime):

    def format(self, value):
        if isinstance(value, datetime):
            return super().format(value.replace(tzinfo=pytz.utc))
        else:
            return str(value)


rfc3339DateTime = RFC3339DateTime()


class DotKeyFieldMixin:
    ''' Allows use of flask_restplus fields with '.' in key names. By default, '.'
    is used as a separator for accessing nested properties. Mixin prevents this,
    allowing fields to use '.' in the key names.

    Example of issue:
    >>> data = {"my.dot.field": 1234}
    >>> model = {"my.dot.field: fields.String}
    >>> marshal(data, model)
    {"my.dot.field:": None}

    flask_restplus tries to fetch values for data['my']['dot']['field'] instead
    of data['my.dot.field'] which is the desired behaviour in this case.
    '''

    def output(self, key, obj, **kwargs):
        transformed_obj = {k.replace(".", "___"): v for k, v in obj.items()}
        transformed_key = key.replace(".", "___")
        # if self.attribute is set and contains '.' super().output() will
        # use '.' as a separator for nested access.
        # -> temporarily set to None to overcome this
        with self.toggle_attribute():
            data = super().output(transformed_key, transformed_obj)
        return data

    @contextmanager
    def toggle_attribute(self):
        ''' Context manager to temporarily set self.attribute to None

        Yields self.attribute before setting to None
        '''
        attribute = self.attribute
        self.attribute = None
        yield attribute
        self.attribute = attribute


class DotKeyNested(DotKeyFieldMixin, fields.Nested):
    pass
