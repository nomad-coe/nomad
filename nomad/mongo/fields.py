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

from datetime import datetime
from typing import TYPE_CHECKING, Any

from mongoengine import DateTimeField

from nomad.utils import normalize_datetime_utc


class UTCDateTimeField(DateTimeField):
    """
    DateTimeField that normalizes values to UTC.

    Read:
    - returns timezone-aware UTC datetime
    """

    if TYPE_CHECKING:

        def __get__(self, instance: Any, owner: type[Any]) -> datetime | None: ...

        def __set__(self, instance: Any, value: datetime | None) -> None: ...

    def to_python(self, value):
        parsed = super().to_python(value)

        if isinstance(parsed, datetime):
            return normalize_datetime_utc(parsed)
        return parsed
