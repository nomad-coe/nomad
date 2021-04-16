from typing import Tuple, List, Union, Dict, Set, Any
from fastapi import HTTPException
from elasticsearch_dsl import Search, Q

from optimade.filterparser import LarkParser
from optimade.server.entry_collections import EntryCollection
from optimade.server.query_params import EntryListingQueryParams, SingleEntryQueryParams
from optimade.server.exceptions import BadRequest
from optimade.server.mappers import StructureMapper
from optimade.models import StructureResource

from nomad import datamodel, files, search, utils
from nomad.normalizing.optimade import (
    optimade_chemical_formula_reduced, optimade_chemical_formula_anonymous,
    optimade_chemical_formula_hill)

from .filterparser import _get_transformer as get_transformer


logger = utils.get_logger(__name__)


class StructureCollection(EntryCollection):

    def __init__(self):
        super().__init__(
            resource_cls=StructureResource,
            resource_mapper=StructureMapper,
            transformer=get_transformer(nomad_properties='dft', without_prefix=False))

        self.parser = LarkParser(version=(1, 0, 0), variant="default")

        # check aliases do not clash with mongo operators
        self._check_aliases(self.resource_mapper.all_aliases())
        self._check_aliases(self.resource_mapper.all_length_aliases())

        self.client = Search

    def _base_search_request(self):
        request = search.SearchRequest().owner('public', None)
        request.search_parameter('processed', True)
        # TODO use the elastic annotations when done
        request.query(Q('exists', field='dft.optimade.elements'))
        return request

    def __len__(self) -> int:
        # TODO cache
        return self._base_search_request().execute()['total']

    def count(self, **kwargs) -> int:
        # This seams solely mongodb specific
        raise NotImplementedError()

    def find(
            self,
            params: Union[EntryListingQueryParams, SingleEntryQueryParams]) \
            -> Tuple[List[StructureResource], int, bool, set]:

        criteria = self.handle_query_params(params)
        single_entry = isinstance(params, SingleEntryQueryParams)
        response_fields = criteria.pop("fields")

        results, data_returned, more_data_available = self._run_db_query(
            criteria, single_entry=isinstance(params, SingleEntryQueryParams)
        )

        results = self._es_to_optimade_results(results, response_fields=response_fields)

        if single_entry:
            results = results[0] if results else None

            if data_returned > 1:
                raise HTTPException(
                    status_code=404,
                    detail=f'Instead of a single entry, {data_returned} entries were found')

        exclude_fields = self.all_fields - response_fields

        return (
            results,
            data_returned,
            more_data_available,
            exclude_fields,
        )

    def _check_aliases(self, aliases):
        pass

    def _es_to_optimade_result(
            self, es_result: dict,
            response_fields: Set[str],
            upload_files_cache: Dict[str, files.UploadFiles]) -> StructureResource:

        calc_id, upload_id = es_result['calc_id'], es_result['upload_id']
        upload_files = upload_files_cache.get(upload_id)

        if upload_files is None:
            upload_files = files.UploadFiles.get(upload_id, is_authorized=lambda: True)
            if upload_files is None:
                logger.error('missing upload', upload_id=upload_id)
                return None

            upload_files_cache[upload_id] = upload_files

        try:
            archive = upload_files.read_archive(calc_id)
        except KeyError:
            logger.error('missing archive entry', upload_id=upload_id, calc_id=calc_id)
            return None

        metadata = archive[calc_id]['section_metadata'].to_dict()
        entry = datamodel.EntryMetadata.m_from_dict(metadata)

        attrs = entry.dft.optimade.m_to_dict()

        attrs['immutable_id'] = calc_id
        attrs['last_modified'] = entry.last_processing if entry.last_processing is not None else entry.upload_time

        # TODO this should be removed, once all data is reprocessed with the right normalization
        attrs['chemical_formula_reduced'] = optimade_chemical_formula_reduced(
            attrs['chemical_formula_reduced'])
        attrs['chemical_formula_anonymous'] = optimade_chemical_formula_anonymous(
            attrs['chemical_formula_reduced'])
        attrs['chemical_formula_hill'] = optimade_chemical_formula_hill(
            attrs['chemical_formula_hill'])
        attrs['chemical_formula_descriptive'] = attrs['chemical_formula_hill']
        dimension_types = attrs['dimension_types']
        if isinstance(dimension_types, int):
            attrs['dimension_types'] = [1] * dimension_types + [0] * (3 - dimension_types)
            attrs['nperiodic_dimensions'] = dimension_types
        elif isinstance(dimension_types, list):
            attrs['nperiodic_dimensions'] = sum(dimension_types)

        if response_fields is not None:
            for request_field in response_fields:
                if not request_field.startswith('_nmd_'):
                    continue

                try:
                    if request_field.startswith('_nmd_dft_'):
                        attrs[request_field] = getattr(entry.dft, request_field[9:])
                    else:
                        attrs[request_field] = getattr(entry, request_field[5:])
                except AttributeError:
                    # if unknown properties where provided, we will ignore them
                    pass

        return self.resource_cls(
            type='structures',
            id=entry.calc_id,
            attributes=attrs,
            relationships=None)

    def _es_to_optimade_results(self, es_results: List[dict], response_fields: Set[str]):
        upload_files_cache: Dict[str, files.UploadFiles] = {}
        optimade_results = []
        try:
            for es_result in es_results:
                optimade_result = self._es_to_optimade_result(
                    es_result, response_fields, upload_files_cache)
                if optimade_result is not None:
                    optimade_results.append(optimade_result)
        finally:
            for upload_files in upload_files_cache.values():
                upload_files.close()

        return optimade_results

    def _run_db_query(self, criteria: Dict[str, Any], single_entry=False):

        sort, order = criteria.get('sort', (('chemical_formula_reduced', 1),))[0]
        sort_quantity = datamodel.OptimadeEntry.m_def.all_quantities.get(sort, None)
        if sort_quantity is None:
            raise BadRequest(detail='Unable to sort on field %s' % sort)
        sort_quantity_a_optimade = sort_quantity.m_get_annotations('optimade')
        if not sort_quantity_a_optimade.sortable:
            raise BadRequest(detail='Unable to sort on field %s' % sort)

        search_request = self._base_search_request().include('calc_id', 'upload_id')

        if criteria.get("filter", False):
            search_request.query(criteria["filter"])

        es_response = search_request.execute_paginated(
            page_offset=criteria.get('skip', 0),
            per_page=criteria['limit'],
            order=order,
            order_by='dft.optimade.%s' % sort)
        results = es_response['results']

        data_returned = es_response['pagination']['total']
        more_data_available = data_returned >= criteria.get('skip', 0) + criteria['limit']

        return results, data_returned, more_data_available

    def insert(self, *args, **kwargs):
        # This is used to insert test records during OPT tests. This should never be necessary
        # on our implementation. We just need to implement it, because its marked as
        # @abstractmethod.
        raise NotImplementedError()
