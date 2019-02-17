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
This module represents calculations in elastic search.
"""

from typing import Iterable, Dict, Tuple, List
from elasticsearch_dsl import Document, InnerDoc, Keyword, Text, Date, \
    Object, Boolean, Search, Integer, Q, A
import elasticsearch.helpers
import ase.data

from nomad import config, datamodel, infrastructure, datamodel, coe_repo, parsing


class AlreadyExists(Exception): pass


class User(InnerDoc):

    @classmethod
    def from_user_popo(cls, user):
        self = cls(user_id=user.id)

        if 'first_name' not in user:
            user = coe_repo.User.from_user_id(user.id).to_popo()

        name = '%s, %s' % (user['last_name'], user['first_name'])
        self.name = name
        self.name_keyword = name

        return self

    user_id = Keyword()
    name = Text()
    name_keyword = Keyword()


class Dataset(InnerDoc):

    @classmethod
    def from_dataset_popo(cls, dataset):
        return cls(
            id=dataset.id,
            doi=dataset.doi.value if dataset.doi is not None else None,
            name=dataset.name)

    id = Keyword()
    doi = Keyword()
    name = Keyword()


class Entry(Document):
    class Index:
        name = config.elastic.index_name

    upload_id = Keyword()
    upload_time = Date()
    calc_id = Keyword()
    calc_hash = Keyword()
    pid = Keyword()
    mainfile = Keyword()
    files = Keyword(multi=True)
    uploader = Object(User)

    with_embargo = Boolean()
    published = Boolean()

    authors = Object(User, multi=True)
    owners = Object(User, multi=True)
    comment = Text()
    references = Keyword()
    datasets = Object(Dataset)

    formula = Keyword()
    atoms = Keyword(multi=True)
    basis_set = Keyword()
    xc_functional = Keyword()
    system = Keyword()
    crystal_system = Keyword()
    spacegroup = Keyword()
    code_name = Keyword()
    code_version = Keyword()

    n_total_energies = Integer()
    n_geometries = Integer()
    geometries = Keyword(multi=True)
    quantities = Keyword(multi=True)

    # def __init__(self, *args, **kwargs):
    #     super().__init__(*args, **kwargs)
    #     self.authors = []
    #     self.owners = []

    @classmethod
    def from_calc_with_metadata(cls, source: datamodel.CalcWithMetadata) -> 'Entry':
        entry = Entry(meta=dict(id=source.calc_id))
        entry.update(source)
        return entry

    def update(self, source: datamodel.CalcWithMetadata) -> None:
        self.upload_id = source.upload_id
        self.upload_time = source.upload_time
        self.calc_id = source.calc_id
        self.calc_hash = source.calc_hash
        self.pid = str(source.pid)
        self.mainfile = source.mainfile
        self.files = source.files
        self.uploader = User.from_user_popo(source.uploader) if source.uploader is not None else None

        self.with_embargo = source.with_embargo
        self.published = source.published
        self.authors = [User.from_user_popo(user) for user in source.coauthors]
        self.owners = [User.from_user_popo(user) for user in source.shared_with]
        if self.uploader is not None:
            self.authors.append(self.uploader)
            self.owners.append(self.uploader)
        self.comment = source.comment
        self.references = [ref.value for ref in source.references]
        self.datasets = [Dataset.from_dataset_popo(ds) for ds in source.datasets]

        self.formula = source.formula
        self.atoms = list(set(source.atoms))
        self.basis_set = source.basis_set
        self.xc_functional = source.xc_functional
        self.system = source.system
        self.crystal_system = source.crystal_system
        self.spacegroup = source.spacegroup
        self.code_name = source.code_name
        self.code_version = source.code_version

        if source.backend is not None:
            quantities = set()
            geometries = set()
            n_total_energies = 0
            n_geometries = 0

            for meta_info, _, value in source.backend._delegate.results.traverse():
                quantities.add(meta_info)
                if meta_info == 'energy_total':
                    n_total_energies += 1
                if meta_info == 'section_system':
                    n_geometries += 1
                if meta_info == 'configuration_raw_gid':
                    geometries.add(value)

            self.geometries = list(geometries)
            self.quantities = list(quantities)
            self.n_total_energies = n_total_energies
            self.n_geometries = n_geometries


def delete_upload(upload_id):
    """ Delete all entries with given ``upload_id`` from the index. """
    index = Entry._default_index()
    Search(index=index).query('match', upload_id=upload_id).delete()


def publish(calcs: Iterable[datamodel.CalcWithMetadata]) -> None:
    """ Update all given calcs with their metadata and set ``publish = True``. """
    def elastic_updates():
        for calc in calcs:
            entry = Entry.from_calc_with_metadata(calc)
            entry.published = True
            yield entry.to_dict(include_meta=True)

    elasticsearch.helpers.bulk(infrastructure.elastic_client, elastic_updates())


default_aggregations = {
    'atoms': len(ase.data.chemical_symbols),
    'system': 10,
    'crystal_system': 10,
    'code_name': len(parsing.parsers),
    'xc_functional': 10,
    'authors': 10
}


def aggregate_search(
        page: int = 0, per_page: int = 10, q: Q = None,
        aggregations: Dict[str, int] = default_aggregations,
        **kwargs) -> Tuple[int, List[dict], Dict[str, Dict[str, int]]]:
    """
    Performs a search and returns paginated search results and aggregation bucket sizes
    based on key quantities.

    Arguments:
        page: The page to return starting with 0
        per_page: Results per page
        q: An *elasticsearch_dsl* query used to further filter the results (via `and`)
        aggregations: A customized list of aggregations to perform. Keys are index fields,
            and values the amount of buckets to return. Only works on *keyword* field.
        **kwargs: Field, value pairs to search for.

    Returns: A tuple with the total hits, an array with the results, an dictionary with
        the aggregation data.
    """

    search = Search()
    if q is not None:
        search.query(q)
    for key, value in kwargs.items():
        if key == 'comment':
            search = search.query(Q('match', **{key: value}))
        else:
            search = search.query(Q('term', **{key: value}))

    for aggregation, size in aggregations.items():
        if aggregation == 'authors':
            search.aggs.bucket(aggregation, A('terms', field='authors.name_keyword', size=size))
        else:
            search.aggs.bucket(aggregation, A('terms', field=aggregation, size=size))

    response = search[page * per_page: (page + 1) * per_page].execute()  # pylint: disable=no-member

    total_results = response.hits.total
    search_results = [hit.to_dict() for hit in response.hits]

    aggregation_results = {
        aggregation: {
            bucket.key: bucket.doc_count
            for bucket in getattr(response.aggregations, aggregation).buckets
        }
        for aggregation in aggregations.keys()
    }

    return total_results, search_results, aggregation_results


def authors(per_page: int = 10, after: str = None, prefix: str = None) -> Tuple[Dict[str, int], str]:
    """
    Returns the name field for all authors with the number of their calculations in
    their natural order.

    The author buckets of :func:`search` is limit to the top 10 author. This function
    in contrast, allows to paginate through all authors.

    Arguments:
        per_page: The number of authors to return
        after: Only return the authors after the given ``name_keyword`` field value
            (for pagination).
        prefix: Used to do a prefix search on authors. Be aware the return authors also
            contain the coauthors of the actual authors with prefix.

    Returns: A tuple with an ordered dict containing author ``name_keyword`` field value
        and number of calculations and the ``name_keyword`` value of the last author
        (to be used as the next ``after`` value).
    """
    composite = dict(
        size=per_page,
        sources=dict(authors=dict(terms=dict(field='authors.name_keyword'))))
    if after is not None:
        composite.update(after=dict(authors=after))

    body = dict(size=0, aggs=dict(authors=dict(composite=composite)))
    if prefix is not None:
        body.update(query=dict(prefix={'authors.name': prefix}))

    response = infrastructure.elastic_client.search(index=config.elastic.index_name, body=body)
    response = response['aggregations']['authors']
    return {
        bucket['key']['authors']: bucket['doc_count']
        for bucket in response['buckets']}, response.get('after_key', {'authors': None})['authors']
