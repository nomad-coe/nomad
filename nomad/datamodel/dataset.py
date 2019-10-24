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

from typing import Dict, Any
import mongoengine as me
from flask_restplus import fields

from nomad.metainfo import MSection, Quantity, Section


class Dataset(MSection):
    """ A Dataset is attached to one or many entries to form a set of data.

    Args:
        dataset_id: The unique identifier for this dataset as a string. It should be
            a randomly generated UUID, similar to other nomad ids.
        name: The human readable name of the dataset as string. The dataset name must be
            unique for the user.
        user_id: The unique user_id of the owner and creator of this dataset. The owner
            must not change after creation.
        doi: The optional Document Object Identifier (DOI) associated with this dataset.
            Nomad can register DOIs that link back to the respective representation of
            the dataset in the nomad UI. This quantity holds the string representation of
            this DOI. There is only one per dataset.
    """
    dataset_id = Quantity(type=str, a_me=dict(primary_key=True))
    name = Quantity(type=str, a_me=dict(index=True))
    user_id = Quantity(type=str, a_me=dict(index=True))
    doi = Quantity(type=str, a_me=dict(index=True))

    @classmethod
    def get(cls, **kwargs):
        dataset = DatasetME.objects(**kwargs).first()
        if dataset is None:
            raise KeyError('Dataset does not exist.')

        dct = dataset.to_mongo().to_dict()
        del(dct['_id'])
        dct['dataset_id'] = dataset.dataset_id
        return Dataset.m_from_dict(dct)


def generate_flask_restplus_model(api, section_def: Section):
    def generate_field(quantity: Quantity):
        field = None
        if quantity.type == int:
            field = fields.Integer
        elif quantity.type == float:
            field = fields.Float
        elif quantity.type == str:
            field = fields.String
        elif quantity.type == bool:
            field = fields.Boolean
        else:
            raise NotImplementedError

        result = field(description=quantity.description)

        if len(quantity.shape) == 0:
            return result
        elif len(quantity.shape) == 1:
            return fields.List(result)
        else:
            raise NotImplementedError

    return api.model(section_def.name, {
        name: generate_field(quantity)
        for name, quantity in section_def.all_quantities.items()
    })


def generate_mongoengine(section_def: Section):
    def generate_field(quantity: Quantity):
        annotation = quantity.m_annotations.get('me', {})
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
        if quantity.m_annotations.get('me', {}).get('index', False)]

    dct: Dict[str, Any] = dict()
    if len(indexes) > 0:
        dct.update(meta=dict(indexes=indexes))
    dct.update(**{
        name: generate_field(quantity)
        for name, quantity in section_def.all_quantities.items()
    })
    return type(section_def.name, (me.Document,), dct)


DatasetME = generate_mongoengine(Dataset.m_def)
