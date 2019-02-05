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
independent from their representation in the different modules
:py:mod:`nomad.processing`, :py:mod:`nomad.coe_repo`, :py:mod:`nomad.files`,
:py:mod:`nomad.search`.

It is not about representing every detail, but those parts that are directly involved in
api, processing, migration, mirroring, or other 'infrastructure' operations.

Transformations between different implementations of the same entity can be build
and used. To ease the number of necessary transformations the classes
:class:`UploadWithMetadata` and :class:`CalcWithMetadata` can act as intermediate
representations. Therefore, implement only transformation from and to these
classes.

To implement a transformation, provide a transformation method in the source
entity class and register it:

.. code-block:: python

    def to_my_target_entity(self):
        target = MyTargetEntity()
        target.property_x = # your transformation code

        return target

    MyTargetEntity.register_mapping(MySourceEntity.to_my_target_entity)

To apply a transformation, use:

.. code-block:: python

    my_target_entity_instance = my_source_entity_instance.to(MyTargetEntity)
"""

from typing import Type, TypeVar, Union, Iterable, cast, Callable, Dict
import datetime

from nomad import utils

T = TypeVar('T')


class Entity():
    """
    A common base class for all nomad entities. It provides the functions necessary
    to apply transformations.
    """
    mappings: Dict[Type['Entity'], Callable[['Entity'], 'Entity']] = dict()

    @classmethod
    def register_mapping(
            cls, from_type: Type['Entity'], mapping: Callable[['Entity'], 'Entity']):
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
        """
        Create an entity of this class from an instance of a different
        class for the same entity type. This requires a registered mapping or
        overloaded implementation.

        Arguments:
            obj: The source entity instance.
        """
        return cls.mappings[obj.__class__](obj)

    def to(self, entity_cls: Type[T]) -> T:
        """
        Either provides a type cast if it already has the right type, or adapt
        the type using the :func:`load_from` of the target class `entity_cls`.
        """
        if (isinstance(self, entity_cls)):
            return cast(T, self)
        else:
            return cast(T, cast(Type[Entity], entity_cls).load_from(self))


class Calc(Entity):

    @property
    def pid(self) -> Union[int, str]:
        """ The repisitory pid, either as int or str. """
        raise NotImplementedError

    @property
    def mainfile(self) -> str:
        """ The mainfile (path) that identifies a calc within an upload. """
        raise NotImplementedError

    @property
    def calc_id(self) -> str:
        """ The internal UUID based on upoad_id and mainfile. """
        raise NotImplementedError

    @property
    def upload(self) -> 'Upload':
        """ A reference to the upload (of same implementation). """
        raise NotImplementedError


class Upload(Entity):

    @property
    def upload_id(self) -> str:
        """ The randomly choosed upload UUID """
        return '<not assigned>'

    @property
    def upload_time(self) -> Type[datetime.datetime]:
        """ The upload time assigned by the API after receiving the uploaded raw archive. """
        raise NotImplementedError

    @property
    def uploader(self):
        """ A reference to the uploaded (i.e. :class:`coe_repo.User`) """
        raise NotImplementedError

    @property
    def calcs(self) -> Iterable[Calc]:
        """ A list of references to the upload's calcs (of same implementation) """
        raise NotImplementedError


class UploadWithMetadata(dict, Entity):
    """
    See :class:`CalcWithMetadata`.
    """

    def __init__(self, upload_id):
        self.upload_id = upload_id


class CalcWithMetadata(utils.POPO, Entity):
    """
    A dict/POPO class that can be used for mapping calc representations with calc metadata.
    We have many representations of calcs and their calc metadata. To avoid implement
    mappings between all combinations, just implement mappings with the class and use
    mapping transitivity. E.g. instead of A -> B, A -> this -> B.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.upload = UploadWithMetadata(kwargs['upload_id'])
