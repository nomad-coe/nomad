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

'''
Adds mongoengine supports to the metainfo. Allows to create, save, and get metainfo
sections from mongoengine. Currently no sub-section support. The annotation key is 'mongo'.
'''

from typing import Any, Dict, List

from .metainfo import DefinitionAnnotation, SectionAnnotation, Annotation, MSection, Datetime, Quantity


class Mongo(DefinitionAnnotation):
    '''
    This annotation class can be used to extend metainfo quantities. It enables and
    details the mapping of quantities to fields in mongoengine documents.

    Attributes:
        index: A boolean indicating that this quantity should be indexed.
        primary_key: A boolean indicating that this quantity is the primary key.
    '''
    def __init__(
            self, index: bool = False, primary_key: bool = False,
            **kwargs):
        self.primary_key = primary_key
        self.index = index
        self.kwargs = kwargs

        if kwargs is None:
            self.kwargs = {}

        if self.primary_key:
            kwargs.update(primary_key=primary_key)


class MongoDocument(SectionAnnotation):
    '''
    This annotation class can be used to extend metainfo section. It allows to get
    the mongoengine document class to store instances of this section in mongodb. It
    also provides access to the respective mongodb collection.
    '''
    def __init__(self):
        self._mongoengine_cls = None

        self.primary_key: Mongo = None

    def new(self, section):
        return dict(mongo=MongoInstance(section))

    @property
    def mongo_cls(self):
        '''
        The mongoengine document class for this section. Only quantities with :class:`Mongo`
        annotation are mapped to fields.
        '''
        if self._mongoengine_cls is not None:
            return self._mongoengine_cls

        import mongoengine as me

        def generate_field(quantity: Quantity, annotation: Mongo):
            field = None
            if quantity.type == int:
                field = me.IntField
            elif quantity.type == float:
                field = me.FloatField
            elif quantity.type == str:
                field = me.StringField
            elif quantity.type == bool:
                field = me.BooleanField
            elif quantity.type == Datetime:
                field = me.DateTimeField
            else:
                raise NotImplementedError

            result = field(default=quantity.default, **annotation.kwargs)

            if len(quantity.shape) == 0:
                return result
            elif len(quantity.shape) == 1:
                return me.ListField(result)
            else:
                raise NotImplementedError

        indexes: List[str] = []
        dct: Dict[str, Any] = {}

        for quantity in self.definition.all_quantities.values():
            annotation = quantity.m_get_annotations(Mongo)
            if annotation is None:
                continue

            if annotation.index:
                indexes.append(quantity.name)

            dct[quantity.name] = generate_field(quantity, annotation)

            if annotation.primary_key:
                self.primary_key = annotation

        if len(indexes) > 0:
            dct['meta'] = dict(indexes=indexes)

        self._mongoengine_cls = type(self.definition.name, (me.Document,), dct)
        return self._mongoengine_cls

    def objects(self, *args, **kwargs):
        '''
        Allows access to the underlying collection objects function.
        Returns mongoengine document instances, not metainfo section instances.
        '''
        return self.mongo_cls.objects(*args, **kwargs)

    def get(self, **kwargs):
        '''
        Returns the first entry that matches the given objects query as metainfo
        section instance. Raises KeyError.
        '''
        mongo_instance = self.objects(**kwargs).first()
        if mongo_instance is None:
            raise KeyError(str(kwargs))
        return self.to_metainfo(mongo_instance)

    def to_metainfo(self, mongo_instance):
        '''
        Turns the given mongoengine document instance into its metainfo section instance
        counterpart.
        '''
        section = self.definition.section_cls()
        section.a_mongo.mongo_instance = mongo_instance
        for name, quantity in self.definition.all_quantities.items():
            if quantity.m_get_annotations(Mongo) is not None:
                value = getattr(mongo_instance, name)
                if value is not None:
                    section.m_set(quantity, value)

        return section


class MongoInstance(Annotation):
    '''
    The annotation that is automatically added to all instances of sections that
    feature the :class:`MongoDocument` annotation.
    '''

    def __init__(self, section: MSection):
        self.section = section
        self.mongo_instance = None

    def save(self):
        ''' Saves the section as mongo entry. Does an upsert. '''
        if self.mongo_instance is None:
            return self.create()

        for quantity_name, quantity in self.section.m_def.all_quantities.items():
            value = self.section.m_get(quantity)

            setattr(self.mongo_instance, quantity_name, value)

        self.mongo_instance.save()
        return self.section

    def create(self):
        ''' Creates a new mongo entry and saves it. '''
        self.mongo_instance = self.section.m_def.a_mongo.mongo_cls()
        return self.save()

    def delete(self):
        ''' Deletes the respective entry from mongodb. '''
        self.mongo_instance.delete()
        self.mongo_instance = None
        return self.section
