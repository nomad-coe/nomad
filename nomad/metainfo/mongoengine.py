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
Adds mongoengine supports to the metainfo. Allows to create, save, and get metainfo
sections from mongoengine. Currently no sub-section support. The annotation key is "a_me",
the annotation object support the following keys:

- ``primary_key``: *Bool*, renders the quantity to be the primary key.
- ``index``: *Bool*, adds this quantity to the index
"""

from typing import Any, Dict
import mongoengine as me

from .metainfo import Section, Quantity, Datetime


def init_section(section_cls):
    section_def = section_cls.m_def
    assert section_def.m_annotations.get('me') is None, 'Can only initialize once'
    section_def.m_annotations['me'] = MESection(section_cls)
    # assert getattr(section_cls, '__init__', None) is None, 'Only section classes without constructor can be used for mongoengine'

    def __init__(self, *args, **kwargs):
        super(section_cls, self).__init__(*args, **kwargs)
        self.m_annotations['me'] = MEInstance(section_cls.m_def.m_x('me'), self)

    section_cls.__init__ = __init__


class MESection():
    def __init__(self, section_cls):
        self.section_cls = section_cls
        self.me_cls = generate_mongoengine(section_cls.m_def)

        section_def = self.section_cls.m_def

        id_quantity = None
        for quantity in section_def.all_quantities.values():
            annotation = quantity.m_annotations.get('me', None)
            if annotation is not None and annotation.get('primary_key', False):
                id_quantity = quantity.name
        assert id_quantity is not None, 'Section %s has no mongoengine primary key' % section_def

        self.id_quantity = id_quantity

    def objects(self, *args, **kwargs):
        return self.me_cls.objects(*args, **kwargs)

    def get(self, **kwargs):
        me_obj = self.objects(**kwargs).first()
        if me_obj is None:
            raise KeyError
        return self.to_metainfo(me_obj)

    def to_metainfo(self, me_obj):
        dct = me_obj.to_mongo().to_dict()
        del(dct['_id'])
        dct[self.id_quantity] = getattr(me_obj, self.id_quantity)
        section = self.section_cls.m_from_dict(dct)  # pylint: disable=no-member
        section.m_x('me').me_obj = me_obj
        return section


class MEInstance():
    def __init__(self, me_section: MESection, metainfo):
        self.me_section = me_section
        self.metainfo = metainfo
        self.me_obj = None

    def save(self):
        if self.me_obj is None:
            return self.create()

        for quantity_name, quantity in self.metainfo.m_def.all_quantities.items():
            me_value = self.metainfo.m_get(quantity)

            setattr(self.me_obj, quantity_name, me_value)

        self.me_obj.save()
        return self.metainfo

    def create(self):
        self.me_obj = self.me_section.me_cls()
        return self.save()

    def delete(self):
        self.me_obj.delete()
        self.me_obj = None
        return self.metainfo


def generate_mongoengine(section_def: Section):
    def generate_field(quantity: Quantity):
        annotation = quantity.m_x('me', {})
        annotation.pop('index', None)

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

        result = field(default=quantity.default, **annotation)

        if len(quantity.shape) == 0:
            return result
        elif len(quantity.shape) == 1:
            return me.ListField(result)
        else:
            raise NotImplementedError

    indexes = [
        quantity.name
        for quantity in section_def.all_quantities.values()
        if quantity.m_annotations.get('a_me', {}).get('index', False)]

    dct: Dict[str, Any] = dict()
    if len(indexes) > 0:
        dct.update(meta=dict(indexes=indexes))
    dct.update(**{
        name: generate_field(quantity)
        for name, quantity in section_def.all_quantities.items()
    })
    return type(section_def.name, (me.Document,), dct)
