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

from typing import Callable, Any, Dict, cast
import uuid
import numpy as np
import pint.quantity


from .metainfo import (
    Section, Quantity, MSection, MEnum, Datetime, Reference, Annotation, SectionAnnotation,
    DefinitionAnnotation)

'''
This module provides metainfo annotation class :class:`Elastic` and
class:`ElasticDocument` that allows to configure metainfo definitions for the use of
metainfo data in elastic search.
'''


class ElasticDocument(SectionAnnotation):
    '''
    This annotation class can be used to extend metainfo sections. It allows to detail
    how section instances (and their sub sections and quantities) should be represented in
    an elastic search index.

    In general sections with this annotation are mapped to elasticsearch_dsl document
    classes, sub sections become inner documents, and quantities with the :class:`Elastic`
    extension become fields in their respective document.

    Attributes:
        document: The elasticsearch_dsl document class that was generated from the
            metainfo section
    '''

    _all_documents: Dict[str, Any] = {}

    def __init__(self, index_name: str = None, id: Callable[[Any], str] = None):
        """
        Args:
            index_name: This is used to optionally add the index_name to the resulting
                elasticsearch_dsl document.
            id: A callable that produces an id from a section instance that is used as id
                for the respective elastic search index entry. The default will be randomly
                generated UUID(4).
        """
        self.index_name = index_name
        self.id = id

        self.fields: Dict[Quantity, str] = {}

    def new(self, section):
        return dict(elastic=ElasticEntry(section))

    @classmethod
    def create_index_entry(cls, section: MSection):
        ''' Creates an elasticsearch_dsl document instance for the given section. '''
        from elasticsearch_dsl import Object

        m_def = section.m_def
        annotation = m_def.m_get_annotations(ElasticDocument)
        document_cls = ElasticDocument._all_documents[m_def.qualified_name()]

        if annotation is None:
            obj = document_cls()
        else:
            if annotation.id is not None:
                id = annotation.id(section)
            else:
                id = uuid.uuid4()
            obj = document_cls(meta=dict(id=id))

        for quantity in m_def.all_quantities.values():
            for annotation in quantity.m_get_annotations(Elastic, as_list=True):
                if annotation.mapping is None:
                    continue

                value = annotation.value(section)
                if value is None or value == []:
                    continue

                # By default the full section is resolved for references
                if isinstance(quantity.type, Reference) and isinstance(annotation.mapping, Object):
                    if quantity.is_scalar:
                        value = ElasticDocument.create_index_entry(cast(MSection, value))
                    else:
                        value = [ElasticDocument.create_index_entry(item) for item in value]

                # Only the magnitude of scalar Pint quantity objects is stored
                if quantity.is_scalar and isinstance(value, pint.quantity._Quantity):
                    value = value.magnitude

                setattr(obj, annotation.field, value)

        for sub_section in m_def.all_sub_sections.values():
            try:
                if sub_section.repeats:
                    mi_values = list(section.m_get_sub_sections(sub_section))
                    if len(mi_values) == 0:
                        continue
                    value = [ElasticDocument.create_index_entry(value) for value in mi_values]
                else:
                    mi_value = section.m_get_sub_section(sub_section, -1)
                    if mi_value is None:
                        continue
                    value = ElasticDocument.create_index_entry(mi_value)

                setattr(obj, sub_section.name, value)
            except KeyError:
                # the sub section definition has no elastic quantities and therefore not
                # corresponding document class
                pass

        return obj

    @classmethod
    def index(cls, section: MSection, **kwargs):
        ''' Adds the given section to its elastic search index. '''
        entry = cls.create_index_entry(section)
        entry.save(**kwargs)
        return entry

    @property
    def document(self):
        return ElasticDocument.create_document(self.definition, index_name=self.index_name)

    @classmethod
    def create_document(
            cls, section: Section, attrs: Dict[str, Any] = None,
            prefix: str = None, index_name: str = None, root=True):
        '''
        Create all elasticsearch_dsl mapping classes for the section and its sub sections.
        '''
        document = cls._all_documents.get(section.qualified_name())
        if document is not None:
            return document

        from elasticsearch_dsl import Document, InnerDoc, Keyword, Date, Integer, Boolean, Object, Double, Float, Long, Nested

        if attrs is None:
            attrs = {}

        # create an field for each sub section
        for sub_section in section.all_sub_sections.values():
            sub_section_prefix = '%s.%s' % (prefix, sub_section.name) if prefix else sub_section.name

            inner_document = ElasticDocument.create_document(
                sub_section.sub_section, prefix=sub_section_prefix, index_name=index_name, root=False)
            if inner_document is not None:
                for annotation in sub_section.m_get_annotations(Elastic, as_list=True):
                    if annotation.nested:
                        assert sub_section.repeats, (
                            "Nested fields should be repeatable. If the subsection cannot be repeated, "
                            "define it as unnested instead."
                        )
                        attrs[sub_section.name] = Nested(inner_document)

                    annotation.register(prefix, annotation.field, index_name)

                if sub_section.name not in attrs:
                    attrs[sub_section.name] = Object(inner_document)

        # create an field for each quantity
        for quantity in section.all_quantities.values():
            first = True
            for annotation in quantity.m_get_annotations(Elastic, as_list=True):
                if annotation.mapping is None and first:
                    kwargs = dict(index=annotation.index)

                    # Find a mapping based on quantity type if not explicitly given
                    if quantity.type == str:
                        annotation.mapping = Keyword(**kwargs)
                    elif quantity.type in [float, np.float64] and quantity.is_scalar:
                        annotation.mapping = Double(**kwargs)
                    elif quantity.type == np.float32 and quantity.is_scalar:
                        annotation.mapping = Float(**kwargs)
                    elif quantity.type in [int, np.int32] and quantity.is_scalar:
                        annotation.mapping = Integer(**kwargs)
                    elif quantity.type == np.int64 and quantity.is_scalar:
                        annotation.mapping = Long(**kwargs)
                    elif quantity.type == bool:
                        annotation.mapping = Boolean(**kwargs)
                    elif quantity.type == Datetime:
                        annotation.mapping = Date(**kwargs)
                    elif isinstance(quantity.type, Reference):
                        inner_prefix = annotation.field
                        if prefix is not None:
                            inner_prefix = '%s.%s' % (prefix, inner_prefix)
                        inner_document = ElasticDocument.create_document(
                            cast(Section, quantity.type.target_section_def), prefix=inner_prefix, index_name=index_name, root=False)
                        annotation.mapping = Object(inner_document)
                    elif isinstance(quantity.type, MEnum):
                        annotation.mapping = Keyword(**kwargs)
                    else:
                        raise NotImplementedError(
                            'Quantity type %s for quantity %s is not supported.' % (quantity.type, quantity))

                assert first or annotation.mapping is None, 'Only the first Elastic annotation is mapped'

                if first:
                    assert annotation.field not in attrs, 'Elastic fields must be unique'
                    attrs[annotation.field] = annotation.mapping
                annotation.register(prefix, annotation.field, index_name)

                first = False

        if len(attrs) == 0:
            # do not create a document/inner document class, if no elastic quantities are defined
            return None

        doc_cls_obj = InnerDoc
        if root:
            doc_cls_obj = Document
        document = type(section.name, (doc_cls_obj,), attrs)
        cls._all_documents[section.qualified_name()] = document
        return document


class ElasticEntry(Annotation):
    def __init__(self, section: MSection):
        self.section = section

    def index(self, **kwargs):
        return ElasticDocument.index(self.section, **kwargs)

    def create_index_entry(self):
        return ElasticDocument.create_index_entry(self.section)


class Elastic(DefinitionAnnotation):
    '''
    This annotation class can be used to extend metainfo quantities. It allows to detail
    how this quantity should be represented in an elastic search index.

    Arguments:
        field:
            The name of the field for this quantity in the elastic search index. The
            default is the quantity name.
        mapping: A valid elasticsearch_dsl mapping. Default is ``Keyword()``.
        value:
            A callable that is applied to the containering section to get a value for
            this quantity when saving the section in the elastic search index. By default
            this will be the serialized quantity value.
        index:
            A boolean that indicates if this quantity should be indexed or merely be
            part of the elastic document ``_source`` without being indexed for search.
        aggregateable:
            A boolean that determines, if this quantity can be used in aggregations
    '''
    def __init__(
            self,
            field: str = None,
            mapping: Any = None,
            value: Callable[[Any], Any] = None,
            index: bool = True):

        self.field = field
        self.mapping = mapping
        self.value = value
        self.index = index

        self.prefix = None
        self.qualified_field = field

    def init_annotation(self, definition):
        super().init_annotation(definition)

        if self.field is None:
            self.field = definition.name

        if self.value is None:
            self.value = lambda section: section.m_get(definition)

    def register(self, prefix: str, field: str, index: str):

        if prefix is None:
            self.qualified_field = field
        else:
            self.qualified_field = '%s.%s' % (prefix, field)

    @property
    def aggregateable(self):
        return self.mapping is None or self.mapping.__class__.__name__ == 'Keyword'
