from typing import Dict, Any
from flask import request, g
from flask_restplus import Resource, fields, abort
import mongoengine as me

from nomad import utils
from nomad.metainfo import MSection, Quantity, Section
from nomad.app.utils import with_logger

from .api import api
from .auth import authenticate
from .common import pagination_model, pagination_request_parser


ns = api.namespace(
    'datasets',
    description='Datasets allow to create sets of related data.')


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


def generate_flask_restplus_model(section_def: Section):
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


dataset_model = generate_flask_restplus_model(Dataset.m_def)
dataset_list_model = api.model('DatasetList', {
    'pagination': fields.Nested(model=pagination_model),
    'results': fields.List(fields.Nested(model=dataset_model, skip_none=True))
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


@ns.route('/')
class DatasetListResource(Resource):
    @api.doc('list_datasets')
    @api.marshal_with(dataset_list_model, skip_none=True, code=200, description='Dateset send')
    @api.expect(pagination_request_parser)
    @authenticate(required=True)
    def get(self):
        """ Retrieve a list of all datasets of the authenticated user. """
        try:
            page = int(request.args.get('page', 1))
            per_page = int(request.args.get('per_page', 10))
        except Exception:
            abort(400, message='bad parameter types')

        result_query = DatasetME.objects(user_id=g.user.user_id)
        return dict(
            pagination=dict(total=result_query.count(), page=page, per_page=per_page),
            results=result_query[(page - 1) * per_page: page * per_page]), 200

    @api.doc('create_dataset')
    @api.response(400, 'The provided data is malformed or a dataset with the name already exists')
    @api.marshal_with(dataset_model, skip_none=True, code=200, description='Dateset send')
    @api.expect(dataset_model)
    @authenticate(required=True)
    def put(self):
        """ Creates a new dataset. """
        data = request.get_json()
        if data is None:
            data = {}

        # unique name
        name = data.get('name', None)
        if name is None:
            abort(400, 'Must provide a dataset name.')

        if DatasetME.objects(user_id=g.user.user_id, name=name).count() > 0:
            abort(400, 'A dataset with name %s does already exist for the current user.' % name)

        # only admin can set user or doi
        if any(key in data for key in ['user_id', 'doi', 'dataset_id']):
            if not g.user.is_admin():
                abort(400, 'The dataset contains information you are not allowed to set.')

        # no other keys
        if any(key not in Dataset.m_def.all_quantities for key in data):
            abort(400, 'The dataset contains unknown keys.')

        if 'user_id' not in data:
            data['user_id'] = g.user.user_id
        dataset_id = data.pop('dataset_id', utils.create_uuid())
        return DatasetME(dataset_id=dataset_id, **data).save(), 200


@ns.route('/<string:name>')
@api.doc(params=dict(name='The name of the requested dataset.'))
class DatasetResource(Resource):
    @api.doc('get_dataset')
    @api.response(404, 'The dataset does not exist')
    @api.marshal_with(dataset_model, skip_none=True, code=200, description='Dateset send')
    @authenticate(required=True)
    def get(self, name: str):
        """ Retrieve a dataset by name. """
        result = DatasetME.objects(user_id=g.user.user_id, name=name).first()
        if result is None:
            abort(404, 'Dataset with name %s does not exist for current user' % name)

        return result

    @api.doc('assign_doi')
    @api.response(404, 'The dataset does not exist')
    @api.response(400, 'The dataset already has a DOI')
    @api.marshal_with(dataset_model, skip_none=True, code=200, description='DOI assigned')
    @authenticate(required=True)
    @with_logger
    def post(self, name: str, logger):
        """ Assign a DOI to the dataset. """
        result = DatasetME.objects(user_id=g.user.user_id, name=name).first()
        if result is None:
            abort(404, 'Dataset with name %s does not exist for current user' % name)

        logger.error('assign datasets is not implemented', user_id=g.user.user_id)

        return result

    @api.doc('delete_dataset')
    @api.response(404, 'The dataset does not exist')
    @api.response(400, 'The dataset has a DOI and cannot be deleted')
    @api.marshal_with(dataset_model, skip_none=True, code=200, description='Dateset deleted')
    @authenticate(required=True)
    def delete(self, name: str):
        """ Assign a DOI to the dataset. """
        result = DatasetME.objects(user_id=g.user.user_id, name=name).first()
        if result is None:
            abort(404, 'Dataset with name %s does not exist for current user' % name)
        if result.doi is not None:
            abort(400, 'Dataset with name %s has a DOI and cannot be deleted' % name)

        result.delete()

        return result
