# Copyright 2020 Markus Scheidgen
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

from rdflib import Graph, Literal, RDF, URIRef, BNode
from rdflib.namespace import Namespace, DCAT, DCTERMS as DCT

from nomad import config
from nomad.datamodel import User

from nomad.datamodel import EntryMetadata, User

from .api import url


VCARD = Namespace('http://www.w3.org/2006/vcard/ns#')


class Mapping():
    def __init__(self):
        self.g = Graph()
        self.g.namespace_manager.bind('dcat', DCAT)
        self.g.namespace_manager.bind('dct', DCT)
        self.g.namespace_manager.bind('vcard', VCARD)

        self.vcards = {}

    def map_entry(self, entry: EntryMetadata):
        dataset = URIRef(url('datasets', entry.calc_id))

        self.g.add((dataset, RDF.type, DCAT.Dataset))
        self.g.add((dataset, DCT.identifier, Literal(entry.calc_id)))
        self.g.add((dataset, DCT.issued, Literal(entry.upload_time)))
        self.g.add((dataset, DCT.modified, Literal(entry.last_processing)))
        self.g.add((dataset, DCAT.landing_page, URIRef('%s/entry/id/%s/%s' % (
            config.gui_url(), entry.upload_id, entry.calc_id))))
        self.g.add((dataset, DCT.title, Literal('unavailable' if entry.formula is None else entry.formula)))
        self.g.add((dataset, DCT.description, Literal('unavailable' if entry.comment is None else entry.comment)))

        self.g.add((dataset, DCT.publisher, self.map_user(entry.uploader)))
        for author in entry.authors:
            self.g.add((dataset, DCT.creator, self.map_user(author)))

        self.g.add((dataset, DCAT.distribution, self.map_distribution(entry, 'api')))
        self.g.add((dataset, DCAT.distribution, self.map_distribution(entry, 'json')))
        self.g.add((dataset, DCAT.distribution, self.map_distribution(entry, 'raw')))

        return dataset

    def map_user(self, user: User):
        # TODO: creator and publisher are FOAF agents, contactPoint is vCard
        vcard = self.vcards.get(user.user_id)
        if vcard is not None:
            return vcard

        user = User.get(user.user_id)
        vcard = BNode()
        self.g.add((vcard, RDF.type, VCARD.Individual))
        self.g.add((vcard, VCARD.givenName, Literal(user.first_name)))
        self.g.add((vcard, VCARD.familyName, Literal(user.last_name)))
        self.g.add((vcard, VCARD.nickName, Literal(user.username)))
        self.g.add((vcard, VCARD.hasEmail, Literal(user.email)))

        self.vcards[user.user_id] = vcard

        return vcard

    def map_distribution(self, entry, dist_kind):
        if dist_kind == 'api':
            # DataService: API
            service = BNode()
            self.g.add((service, RDF.type, DCAT.DataService))
            self.g.add((service, DCT.title, Literal('NOMAD API')))  # How to include terms from swagger document here?
            self.g.add((service, DCT.description, Literal('Official NOMAD API')))  # same question
            self.g.add((service, DCAT.endpointURL, URIRef('https://nomad-lab.eu/prod/rae/api/')))  # config.api_url() ?
            # not sure if the following needs to be dataset specific:
            self.g.add((service, DCAT.endpointDescription, URIRef('https://nomad-lab.eu/prod/rae/api/swagger.json')))

            # Distribution over API
            dist = BNode()
            self.g.add((dist, DCT.title, Literal('unavailable' if entry.formula is None else entry.formula + 'api')))
            self.g.add((dist, RDF.type, DCAT.Distribution))
            self.g.add((dist, DCAT.accessService, service))
        elif dist_kind == 'json':
            # Distribution as JSON
            dist = BNode()
            self.g.add((dist, RDF.type, DCAT.Distribution))
            self.g.add((dist, DCT.title, Literal('unavailable' if entry.formula is None else entry.formula + 'json')))
            self.g.add((dist, DCAT.mediaType, URIRef('https://www.iana.org/assignments/media-types/application/json')))
            self.g.add((dist, DCAT.packageFormat, URIRef('https://www.iana.org/assignments/media-types/application/zip')))
            self.g.add((dist, DCAT.downloadURL, URIRef(
                'http://nomad-lab.eu/prod/rae/api/archive/download?upload_id=%s&calc_id=%s' % (entry.upload_id, entry.calc_id))))
            self.g.add((dist, DCAT.accessURL, URIRef('%s/entry/id/%s/%s' % (
                config.gui_url(), entry.upload_id, entry.calc_id))))
        elif dist_kind == 'raw':
            # Distribution of the raw data
            dist = BNode()
            self.g.add((dist, RDF.type, DCAT.Distribution))
            self.g.add((dist, DCT.title, Literal('unavailable' if entry.formula is None else entry.formula + 'raw')))
            self.g.add((dist, DCAT.accessURL, URIRef('https://nomad-lab.eu/prod/rae/api/raw/calc/%s/%s' % (
                entry.upload_id, entry.calc_id))))
            self.g.add((dist, DCAT.packageFormat, URIRef('https://www.iana.org/assignments/media-types/application/zip')))

        return dist
