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

from typing import Type, TypeVar, Union, Iterable, cast, Callable, Dict
import datetime

T = TypeVar('T')


class Entity():
    @classmethod
    def load_from(cls: Type[T], obj) -> T:
        raise NotImplementedError

    def to(self, entity_cls: Type[T]) -> T:
        """
        Either provides a type cast if it already has the right type, or adapt
        the type using the :func:`load_from` of the target class :param:`entity_cls`.
        """
        if (isinstance(self, entity_cls)):
            return cast(T, self)
        else:
            return cast(T, cast(Type[Entity], entity_cls).load_from(self))


class Calc(Entity):
    """
    A nomad calculation.

    Attributes:
        pid: The persistent id (pid) for the calculation
        mainfile: The mainfile path relative to upload root
        calc_id: A unique id/checksum that describes unique calculations
        upload: The upload object that this calculation belongs to.
    """
    @property
    def pid(self) -> Union[int, str]:
        raise NotImplementedError

    @property
    def mainfile(self) -> str:
        raise NotImplementedError

    @property
    def calc_id(self) -> str:
        raise NotImplementedError

    @property
    def upload(self) -> 'Upload':
        raise NotImplementedError


class Upload(Entity):
    """
    A nomad upload.

    Attributes:
        upload_id(str): The unique random id that each upload has
        upload_time(datatime): The upload time
        uploader(repo.User): The user that uploaded this upload
        calcs(Iterable[Calc]): An iterable over the calculations of this upload
    """
    @property
    def upload_id(self) -> str:
        return '<not assigned>'

    @property
    def upload_time(self) -> Type[datetime.datetime]:
        raise NotImplementedError

    @property
    def uploader(self):
        raise NotImplementedError

    @property
    def calcs(self) -> Iterable[Calc]:
        raise NotImplementedError


class UploadWithMetadata(dict, Entity):

    def __init__(self, upload_id):
        self.upload_id = upload_id


class CalcWithMetadata(dict, Entity):
    """
    A dict/POPO class that can be used for mapping calc representations with calc metadata.
    We have many representations of calcs and their calc metadata. To avoid implement
    mappings between all combinations, just implement mappings with the class and use
    mapping transitivity. E.g. instead of A -> B, A -> this -> B.

    The other calc representations can register mappings from them, in order to allow
    to use this classes `load_from` method.
    """
    mappings: Dict[Type[Entity], Callable[[Entity], 'CalcWithMetadata']] = dict()

    @classmethod
    def register_mapping(
            cls, from_type: Type[T], mapping: Callable[[T], 'CalcWithMetadata']):
        """
        Register a mapping from instances of another calc representation to instances of
        :class:`CalcWithMetadata`.
        Arguments:
            from_type: The source calc type of the mapping.
            mapping: The mapping itself as a callable that takes a source object of the
                source calc type and returns an instance of :class:`CalcWithMetadata`.
        """
        cls.mappings[from_type] = mapping

    @classmethod
    def load_from(cls, obj):
        return CalcWithMetadata.mappings[obj.__class__](obj)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.upload = UploadWithMetadata(kwargs['upload_id'])

    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: " + name)

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError("No such attribute: " + name)
