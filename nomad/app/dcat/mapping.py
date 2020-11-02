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

        return dataset

    def map_user(self, user: User):
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
