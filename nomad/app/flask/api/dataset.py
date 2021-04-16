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

from flask import request, g
from flask_restplus import Resource, fields, abort
import re

from nomad import utils, processing as proc
from nomad.datamodel import Dataset
from nomad.metainfo.flask_extension import generate_flask_restplus_model
from nomad.doi import DOI

from .. import common
from .api import api
from .auth import authenticate
from .common import pagination_model, pagination_request_parser
from .repo import edit


ns = api.namespace(
    'datasets',
    description='Datasets allow to create sets of related data.')

dataset_model = generate_flask_restplus_model(api, Dataset.m_def)
dataset_list_model = api.model('DatasetList', {
    'pagination': fields.Nested(model=pagination_model, skip_none=True),
    'results': fields.List(fields.Nested(model=dataset_model, skip_none=True))
})

list_datasets_parser = pagination_request_parser.copy()
list_datasets_parser.add_argument('prefix', help='Only return dataset with names that start with prefix.')


@ns.route('/')
class DatasetListResource(Resource):
    @api.doc('list_datasets')
    @api.marshal_with(dataset_list_model, skip_none=True, code=200, description='Dateset send')
    @api.expect(list_datasets_parser)
    @authenticate(required=True)
    def get(self):
        ''' Retrieve a list of all datasets of the authenticated user. '''
        args = {
            key: value for key, value in list_datasets_parser.parse_args().items()
            if value is not None}

        page = args.get('page', 1)
        per_page = args.get('per_page', 10)
        prefix = args.get('prefix', '')

        query_params = dict(user_id=g.user.user_id)
        if prefix != '':
            query_params.update(name=re.compile('^%s.*' % prefix, re.IGNORECASE))

        result_query = Dataset.m_def.a_mongo.objects(**query_params)

        return dict(
            pagination=dict(total=result_query.count(), page=page, per_page=per_page),
            results=[
                Dataset.m_def.a_mongo.to_metainfo(result)
                for result in result_query[(page - 1) * per_page: page * per_page]]), 200

    @api.doc('create_dataset')
    @api.response(400, 'The provided data is malformed or a dataset with the name already exists')
    @api.marshal_with(dataset_model, skip_none=True, code=200, description='Dateset send')
    @api.expect(dataset_model)
    @authenticate(required=True)
    def put(self):
        ''' Creates a new dataset. '''
        data = request.get_json()
        if data is None:
            data = {}

        # unique name
        name = data.get('name', None)
        if name is None:
            abort(400, 'Must provide a dataset name.')

        if Dataset.m_def.a_mongo.objects(user_id=g.user.user_id, name=name).count() > 0:
            abort(400, 'A dataset with name %s does already exist for the current user.' % name)

        # only admin can set user or doi
        if any(key in data for key in ['user_id', 'doi', 'dataset_id']):
            if not g.user.is_admin():
                abort(400, 'The dataset contains information you are not allowed to set.')

        # no other keys
        if any(key not in Dataset.m_def.all_quantities for key in data):  # pylint: disable=all
            abort(400, 'The dataset contains unknown keys.')

        if 'user_id' not in data:
            data['user_id'] = g.user.user_id
        dataset_id = data.pop('dataset_id', utils.create_uuid())
        return Dataset(dataset_id=dataset_id, **data).a_mongo.create(), 200


@ns.route('/<path:name>')
@api.doc(params=dict(name='The name of the requested dataset.'))
class DatasetResource(Resource):
    @api.doc('get_dataset')
    @api.response(404, 'The dataset does not exist')
    @api.marshal_with(dataset_model, skip_none=True, code=200, description='Dateset send')
    @authenticate(required=True)
    def get(self, name: str):
        ''' Retrieve a dataset by name. '''
        try:
            result = Dataset.m_def.a_mongo.get(user_id=g.user.user_id, name=name)
        except KeyError:
            abort(404, 'Dataset with name %s does not exist for current user' % name)

        return result

    @api.doc('assign_doi')
    @api.response(404, 'The dataset does not exist')
    @api.response(400, 'The dataset already has a DOI')
    @api.marshal_with(dataset_model, skip_none=True, code=200, description='DOI assigned')
    @authenticate(required=True)
    def post(self, name: str):
        ''' Assign a DOI to the dataset. '''
        try:
            result = Dataset.m_def.a_mongo.get(user_id=g.user.user_id, name=name)
        except KeyError:
            abort(404, 'Dataset with name %s does not exist for current user' % name)

        if result.doi is not None:
            abort(400, 'Dataset with name %s already has a DOI' % name)

        # check if the DOI can be created
        published_values = proc.Calc.objects(metadata__datasets=result.dataset_id).distinct('metadata.published')

        if False in published_values:
            abort(400, 'Dataset must not contain non published entries.')

        if True not in published_values:
            abort(400, 'Dataset must not be empty.')

        # set the DOI
        doi = DOI.create(title='NOMAD dataset: %s' % result.name, user=g.user)
        doi.create_draft()
        doi.make_findable()

        result.doi = doi.doi

        result.a_mongo.save()
        if doi.state != 'findable':
            common.logger.warning(
                'doi was created, but is not findable', doi=doi.doi, doi_state=doi.state,
                dataset=result.dataset_id)

        # update all affected calcs in the search index
        edit(dict(dataset_id=result.dataset_id))

        return result

    @api.doc('delete_dataset')
    @api.response(404, 'The dataset does not exist')
    @api.response(400, 'The dataset has a DOI and cannot be deleted')
    @api.marshal_with(dataset_model, skip_none=True, code=200, description='Dateset deleted')
    @authenticate(required=True)
    def delete(self, name: str):
        ''' Delete the dataset. '''
        try:
            result = Dataset.m_def.a_mongo.get(user_id=g.user.user_id, name=name)
        except KeyError:
            abort(404, 'Dataset with name %s does not exist for current user' % name)

        if result.doi is not None and len(result.doi) > 0:
            if g.user.is_admin:
                DOI.objects(doi__in=result.doi).delete()
            else:
                abort(400, 'Dataset with name %s has a DOI and cannot be deleted' % name)

        # edit all affected entries
        edit(
            dict(dataset_id=result.dataset_id),
            {'__raw__': {'$pull': {'metadata.datasets': result.dataset_id}}})

        # delete the dataset
        result.a_mongo.delete()

        return result


@ns.route('/doi/<path:doi>')
class RepoPidResource(Resource):
    @api.doc('resolve_doi')
    @api.response(404, 'Dataset with DOI does not exist')
    @api.marshal_with(dataset_model, skip_none=True, code=200, description='DOI resolved')
    @authenticate()
    def get(self, doi: str):
        dataset_me = Dataset.m_def.a_mongo.objects(doi=doi).first()
        if dataset_me is None:
            abort(404, 'Dataset with DOI %s does not exist' % doi)

        return dataset_me, 200
