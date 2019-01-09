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

"""
This module contains classes that allow to represent the core
nomad data entities :class:`Upload` and :class:`Calc` on a high level of abstraction
independent from their representation in the different modules :py:mod:`nomad.repo`,
:py:mod:`nomad.processing`, :py:mod:`nomad.coe_repo`, :py:mod:`nomad.files`.
It is not about representing every detail, but those parts that are directly involved in
api, processing, migration, mirroring, or other 'infrastructure' operations.
"""

from typing import Type, TypeVar, Union, Iterable, cast
import datetime

T = TypeVar('T')


class Entity():
    @classmethod
    def create_from(cls: Type[T], obj) -> T:
        raise NotImplementedError

    def to(self, entity_cls: Type[T]) -> T:
        """
        Either provides a type cast if it already has the right type, or adapt
        the type using the :func:`create_from` of the target class :param:`entity_cls`.
        """
        if (isinstance(self, entity_cls)):
            return cast(T, self)
        else:
            return cast(T, cast(Type[Entity], entity_cls).create_from(self))


class Calc(Entity):
    """
    A nomad calculation.

    Attributes:
        pid: The persistent id (pid) for the calculation
        mainfile: The mainfile path relative to upload root
        calc_hash: A unique hash/checksum that describes unique calculations
        upload: The upload object that this calculation belongs to.
    """
    @property
    def pid(self) -> Union[int, str]:
        raise NotImplementedError

    @property
    def mainfile(self) -> str:
        raise NotImplementedError

    @property
    def calc_hash(self) -> str:
        raise NotImplementedError

    @property
    def upload(self) -> 'Upload':
        raise NotImplementedError


class Upload(Entity):
    """
    A nomad upload.

    Attributes:
        upload_id(str): The unique random id that each upload has
        upload_hash(str): The hash/checksum that describes unique uploads
        upload_time(datatime): The upload time
        uploader(repo.User): The user that uploaded this upload
        calcs(Iterable[Calc]): An iterable over the calculations of this upload
    """
    @property
    def upload_id(self) -> str:
        return '<not assigned>'

    @property
    def upload_hash(self) -> str:
        raise NotImplementedError

    @property
    def upload_time(self) -> Type[datetime.datetime]:
        raise NotImplementedError

    @property
    def uploader(self):
        raise NotImplementedError

    @property
    def calcs(self) -> Iterable[Calc]:
        raise NotImplementedError
