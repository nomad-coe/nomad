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
The repository API of the nomad@FAIRDI APIs. Currently allows to resolve repository
meta-data as well as raw-calculation files.
"""

import os.path
from contextlib import contextmanager
from zipfile import ZIP_DEFLATED

import zipstream
from elasticsearch.exceptions import NotFoundError
from flask import Response, g, request, send_file
from flask_restful import Resource, abort
from werkzeug.exceptions import HTTPException

from nomad.files import RepositoryFile
from nomad.repo import RepoCalc
from nomad.utils import get_logger

from .app import api, app, auth, base_path


class RepoCalcRes(Resource):
    def get(self, upload_hash, calc_hash):
        """
        Get calculation data in repository form, which only entails the quanties shown
        in the repository. This is basically the elastic search index entry for the
        requested calculations. Calcs are references via *upload_hash*, *calc_hash*
        pairs.

        .. :quickref: repo; Get calculation data in repository form.

        **Example request**:

        .. sourcecode:: http

            GET /nomad/api/repo/W36aqCzAKxOCfIiMFsBJh3nHPb4a/7ddvtfRfZAvc3Crr7jOJ8UH0T34I HTTP/1.1
            Accept: application/json

        **Example response**:

        .. sourcecode:: http

            HTTP/1.1 200 OK
            Vary: Accept
            Content-Type: application/json

            {
                "calc_hash":"7ddvtfRfZAvc3Crr7jOJ8UH0T34I",
                "upload_time":"2018-08-30T08:41:51.771367",
                "upload_id":"5b87adb813a441000a70a968",
                "upload_hash":"W36aqCzAKxOCfIiMFsBJh3nHPb4a",
                "mainfile":"RopD3Mo8oMV_-E5bh8uW5PiiCRkH1/data/BrK_svSi/TFCC010.CAB/vasprun.xml.relax1",
                "program_name":"VASP",
                "program_version":"4.6.35  3Apr08 complex  parallel LinuxIFC",
                "chemical_composition_bulk_reduced":"BrKSi2",
                "program_basis_set_type":"plane waves",
                "atom_species":[
                    35,
                    19,
                    14,
                    14
                ],
                "system_type":"Bulk",
                "crystal_system":"orthorhombic",
                "space_group_number":47,
                "configuration_raw_gid":"sq6wTJjRKb2VTajoDLVWDxHCgyN6i",
                "XC_functional_name":"GGA_X_PBE"
            }

        :param string upload_hash: the hash of the upload (from uploaded file contents)
        :param string calc_hash: the hash of the calculation (from mainfile)
        :resheader Content-Type: application/json
        :status 200: calc successfully retrieved
        :status 404: calc with given hashes does not exist
        :returns: the repository calculation entry
        """
        try:
            return RepoCalc.get(id='%s/%s' % (upload_hash, calc_hash)).json_dict, 200
        except NotFoundError:
            abort(404, message='There is no calculation for %s/%s' % (upload_hash, calc_hash))
        except Exception as e:
            abort(500, message=str(e))


class RepoCalcsRes(Resource):
    @auth.login_required
    def get(self):
        """
        Get *'all'* calculations in repository from, paginated.

        .. :quickref: repo; Get *'all'* calculations in repository from, paginated.

        **Example request**:

        .. sourcecode:: http

            GET /nomad/api/repo?page=1&per_page=25 HTTP/1.1
            Accept: application/json

        **Example response**:

        .. sourcecode:: http

            HTTP/1.1 200 OK
            Vary: Accept
            Content-Type: application/json

            {
                "pagination":{
                    "total":1,
                    "page":1,
                    "per_page":25
                },
                "results":[
                    {
                        "calc_hash":"7ddvtfRfZAvc3Crr7jOJ8UH0T34I",
                        "upload_time":"2018-08-30T08:41:51.771367",
                        "upload_id":"5b87adb813a441000a70a968",
                        "upload_hash":"W36aqCzAKxOCfIiMFsBJh3nHPb4a",
                        "mainfile":"RopD3Mo8oMV_-E5bh8uW5PiiCRkH1/data/BrK_svSi/TFCC010.CAB/vasprun.xml.relax1",
                        "program_name":"VASP",
                        "program_version":"4.6.35  3Apr08 complex  parallel LinuxIFC",
                        "chemical_composition_bulk_reduced":"BrKSi2",
                        "program_basis_set_type":"plane waves",
                        "atom_species":[
                            35,
                            19,
                            14,
                            14
                        ],
                        "system_type":"Bulk",
                        "crystal_system":"orthorhombic",
                        "space_group_number":47,
                        "configuration_raw_gid":"sq6wTJjRKb2VTajoDLVWDxHCgyN6i",
                        "XC_functional_name":"GGA_X_PBE"
                    }
                ]
            }

        :qparam int page: the page starting with 1
        :qparam int per_page: desired calcs per page
        :qparam string owner: specifies which cals to return: all|user|staging, default is all
        :resheader Content-Type: application/json
        :status 200: calcs successfully retrieved
        :returns: a list of repository entries in ``results`` and pagination info
        """
        # TODO use argparse? bad request reponse an bad params, pagination as decorator
        page = int(request.args.get('page', 1))
        per_page = int(request.args.get('per_page', 10))
        owner = request.args.get('owner', 'all')

        try:
            assert page >= 1
            assert per_page > 0
        except AssertionError:
            abort(400, message='invalid pagination')

        if owner == 'all':
            search = RepoCalc.search().query('match_all')
        elif owner == 'user':
            if g.user is None:
                abort(401, message='Authentication required for owner value user.')
            search = RepoCalc.search().query('match_all')
            search = search.filter('term', user_id=g.user.email)
        elif owner == 'staging':
            if g.user is None:
                abort(401, message='Authentication required for owner value user.')
            search = RepoCalc.search().query('match_all')
            search = search.filter('term', user_id=g.user.email).filter('term', staging=True)
        else:
            abort(400, message='Invalid owner value. Valid values are all|user|staging, default is all')

        search = search[(page - 1) * per_page: page * per_page]
        return {
            'pagination': {
                'total': search.count(),
                'page': page,
                'per_page': per_page
            },
            'results': [result.json_dict for result in search]
        }


@app.route('%s/raw/<string:upload_hash>/<string:calc_hash>' % base_path, methods=['GET'])
def get_raw(upload_hash, calc_hash):
    """
    Get calculation mainfile raw data. Calcs are references via *upload_hash*, *calc_hash*
    pairs. Returns the mainfile, unless an aux_file is specified. Aux files are stored
    in repository entries. See ``/repo`` endpoint.

    .. :quickref: repo; Get calculation raw data.

    **Example request**:

    .. sourcecode:: http

        GET /nomad/api/raw/W36aqCzAKxOCfIiMFsBJh3nHPb4a/7ddvtfRfZAvc3Crr7jOJ8UH0T34I HTTP/1.1
        Accept: application/gz

    :param string upload_hash: the hash of the upload (from uploaded file contents)
    :param string calc_hash: the hash of the calculation (from mainfile)
    :qparam str auxfile: an optional aux_file to download the respective aux file, default is mainfile
    :qparam all: set any value to get a .zip with main and aux files instead of an individual file
    :resheader Content-Type: application/json
    :status 200: calc raw data successfully retrieved
    :status 404: calc with given hashes does not exist or the given aux file does not exist
    :returns: the raw data in body
    """
    archive_id = '%s/%s' % (upload_hash, calc_hash)
    logger = get_logger(__name__, endpoint='raw', action='get', archive_id=archive_id)

    try:
        repo = RepoCalc.get(id=archive_id)
    except NotFoundError:
        abort(404, message='There is no calculation for %s/%s' % (upload_hash, calc_hash))
    except Exception as e:
        abort(500, message=str(e))

    repository_file = RepositoryFile(upload_hash)

    @contextmanager
    def raw_file(filename):
        try:
            the_file = repository_file.get_file(filename)
            with the_file.open() as f:
                yield f
        except KeyError:
            abort(404, message='The file %s does not exist.' % filename)
        except FileNotFoundError:
            abort(404, message='The file %s does not exist.' % filename)

    get_all = request.args.get('all', None) is not None
    if get_all:
        # retrieve the 'whole' calculation, meaning the mainfile and all aux files as
        # a .zip archive
        def generator():
            """ Stream a zip file with all files using zipstream. """
            def iterator():
                """ Replace the directory based iter of zipstream with an iter over all raw files. """
                def write(filename):
                    """ Write a raw file to the zipstream. """
                    def iter_content():
                        """ Iterate the raw file contents. """
                        with raw_file(filename) as file_object:
                            while True:
                                data = file_object.read(1024)
                                if not data:
                                    break
                                yield data
                    return dict(arcname=filename, iterable=iter_content())

                yield write(repo.mainfile)
                try:
                    for auxfile in repo.aux_files:
                        yield write(os.path.join(os.path.dirname(repo.mainfile), auxfile))
                except Exception as e:
                    logger.error('Exception while accessing auxfiles.', exc_info=e)

            zip_stream = zipstream.ZipFile(mode='w', compression=ZIP_DEFLATED)
            zip_stream.paths_to_write = iterator()

            for chunk in zip_stream:
                yield chunk

        response = Response(generator(), mimetype='application/zip')
        response.headers['Content-Disposition'] = 'attachment; filename={}'.format('%s.zip' % archive_id)
        return response
    else:
        # retrieve an individual raw file
        auxfile = request.args.get('auxfile', None)
        if auxfile:
            filename = os.path.join(os.path.dirname(repo.mainfile), auxfile)
        else:
            filename = repo.mainfile

        try:
            with raw_file(filename) as f:
                rv = send_file(
                    f,
                    mimetype='application/octet-stream',
                    as_attachment=True,
                    attachment_filename=os.path.basename(filename))
                return rv
        except HTTPException as e:
            raise e
        except Exception as e:
            logger = get_logger(
                __name__, endpoint='archive', action='get',
                upload_hash=upload_hash, calc_hash=calc_hash)
            logger.error('Exception on accessing archive', exc_info=e)
            abort(500, message='Could not accessing the archive.')


api.add_resource(RepoCalcsRes, '%s/repo' % base_path)
api.add_resource(RepoCalcRes, '%s/repo/<string:upload_hash>/<string:calc_hash>' % base_path)
