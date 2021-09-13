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
sections from mongoengine. The annotation key is 'mongo'.
'''

from typing import Any, Dict, List

from .metainfo import (
    DefinitionAnnotation, SectionAnnotation, Annotation, MSection, Datetime, Quantity,
    MEnum, JSON)


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
        self.primary_key_name: str = None

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
            elif isinstance(quantity.type, MEnum):
                field = me.StringField
            elif quantity.type == JSON:
                field = me.DictField
            else:
                raise NotImplementedError

            result = field(default=quantity.default, **annotation.kwargs)

            if len(quantity.shape) == 0:
                return result
            elif len(quantity.shape) == 1:
                return me.ListField(result, default=None)
            else:
                raise NotImplementedError

        def create_model_recursive(section, level):
            indexes: List[str] = []
            dct: Dict[str, Any] = {}

            # Add quantities to model
            for quantity in section.all_quantities.values():
                annotation = quantity.m_get_annotations(Mongo)
                if annotation is None:
                    continue

                if annotation.index:
                    indexes.append(quantity.name)

                # Primary key is only stored from the root document.
                if level == 0:
                    if annotation.primary_key:
                        self.primary_key = annotation
                        self.primary_key_name = quantity.name

                dct[quantity.name] = generate_field(quantity, annotation)

            # Add subsections to the model
            for subsection in section.all_sub_sections.values():

                annotation = subsection.sub_section.m_get_annotations(MongoDocument)
                if annotation is None:
                    continue

                embedded_doc_field = type(subsection.sub_section.name, (me.EmbeddedDocumentField,), {})
                model = create_model_recursive(subsection.sub_section, level + 1)
                if subsection.repeats:
                    dct[subsection.name] = me.ListField(embedded_doc_field(model))
                else:
                    dct[subsection.name] = embedded_doc_field(model)

            # Add meta dictionary. The strict mode is set to false in order to
            # not raise an exception when reading data that is not specified in
            # the model.
            meta = {
                "strict": False
            }
            if len(indexes) > 0:
                meta["indexes"] = indexes
            dct['meta'] = meta

            # Return final model
            if level == 0:
                model = type(section.name, (me.Document,), dct)
            else:
                model = type(section.name, (me.EmbeddedDocument,), dct)

            return model

        self._mongoengine_cls = create_model_recursive(self.definition, 0)
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
        section_cls = self.definition.section_cls

        # Get the mongo instance data as dict. This is easy to de-serialize
        # into a metainfo section. If a primary key has been declared, rename
        # the _id field.
        mongo_dict = mongo_instance.to_mongo().to_dict()
        if self.primary_key_name is not None:
            mongo_dict[self.primary_key_name] = mongo_dict["_id"]
        del mongo_dict["_id"]

        section = section_cls.m_from_dict(mongo_dict)
        section.a_mongo.mongo_instance = mongo_instance

        return section


class MongoInstance(Annotation):
    '''
    The annotation that is automatically added to all instances of sections that
    feature the :class:`MongoDocument` annotation.
    '''
    def __init__(self, section: MSection):
        self.section = section
        self.mongo_instance = None
        self._id = None

    def save(self):
        ''' Saves the section as mongo entry. Does an upsert. '''

        # The best way to update a complex entry with mongoengine is to create
        # a new Document instance and specify the target ID which should be
        # updated. The targen ID is taken from an old previously saved
        # instance. If no previous saves have been done, a new object will be
        # created. See discussion at:
        # https://stackoverflow.com/questions/19002469/update-a-mongoengine-document-using-a-python-dict
        data = self.section.m_to_dict()
        if self.mongo_instance is not None:
            data["id"] = self.mongo_instance.id
        mongo_instance = self.section.m_def.a_mongo.mongo_cls(**data, _created=False)
        self.mongo_instance = mongo_instance.save()

        return self.section

    def create(self):
        ''' Creates a new mongo entry and saves it. '''
        return self.save()

    def delete(self):
        ''' Deletes the respective entry from mongodb. '''
        self.mongo_instance.delete()
        self.mongo_instance = None
        return self.section
