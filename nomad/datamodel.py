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
This module contains classes that allow to represent the core
nomad data entities :class:`Upload` and :class:`Calc` on a high level of abstraction
independent from their representation in the different modules
:py:mod:`nomad.processing`, :py:mod:`nomad.coe_repo`, :py:mod:`nomad.parsing`,
:py:mod:`nomad.search`, :py:mod:`nomad.api`, :py:mod:`nomad.migration`.

It is not about representing every detail, but those parts that are directly involved in
api, processing, migration, mirroring, or other 'infrastructure' operations.

Transformations between different implementations of the same entity can be build
and used. To ease the number of necessary transformations the classes
:class:`UploadWithMetadata` and :class:`CalcWithMetadata` can act as intermediate
representations. Therefore, implement only transformation from and to these
classes. These are the implemented transformations:

.. image:: datamodel_transformations.png
"""

from typing import Iterable, List
import datetime

from nomad import utils


class UploadWithMetadata():
    """
    See :class:`CalcWithMetadata`.
    """

    def __init__(self, **kwargs):
        self.upload_id: str = None
        self.uploader: utils.POPO = None
        self.upload_time: datetime.datetime = None

        self.calcs: Iterable[CalcWithMetadata] = list()

        for key, value in kwargs.items():
            setattr(self, key, value)


class CalcWithMetadata():
    """
    A dict/POPO class that can be used for mapping calc representations with calc metadata.
    We have many representations of calcs and their calc metadata. To avoid implement
    mappings between all combinations, just implement mappings with the class and use
    mapping transitivity. E.g. instead of A -> B, A -> this -> B.

    Attributes:
        upload_id: The ``upload_id`` of the calculations upload (random UUID).
        calc_id: The unique mainfile based calculation id.
        upload_time: The time when the calc was uploaded.
        calc_hash: The raw file content based checksum/hash of this calculation.
        pid: The unique persistent id of this calculation.
        mainfile: The upload relative mainfile path.
        files: A list of all files, relative to upload.
        uploader: An object describing the uploading user, has at least ``user_id``
        with_embargo: Show if user set an embargo on the calculation.
        coauthors: List of coauther user objects with at ``user_id``.
        shared_with: List of users this calcs ownership is shared with, objects with at ``user_id``.
        comment: String comment.
        references: Objects describing user provided references, keys are ``id`` and ``value``.
        datasets: Objects describing the datasets, keys are ``id``, ``name``, ``doi``.
            DOI is optional, is an object with key ``id``, ``value``.
        formula: The chemical formula
        atoms: A list of all atoms, as labels. All atoms means the whole composition, with atom labels repeated.
        basis_set: The basis set type of this calculation.
        xc_functional: The class of functional used.
        system: The system type, e.g. Atom/Molecule, 2D, Bulk(3D)
        crystal_system: The symmetry describing crystal_system type.
        spacegroup: The spacegroup, as spacegroup number.
        code_name: The name of the used code.
        code_version: The version of the used code.
    """
    def __init__(self, **kwargs):
        self.upload_id: str = None
        self.calc_id: str = None

        self.upload_time: datetime.datetime = None
        self.calc_hash: str = None
        self.pid: int = None
        self.mainfile: str = None
        self.files: List[str] = None
        self.uploader: utils.POPO = None

        self.with_embargo: bool = None
        self.coauthors: List[utils.POPO] = []
        self.shared_with: List[utils.POPO] = []
        self.comment: str = None
        self.references: List[utils.POPO] = []
        self.datasets: List[utils.POPO] = []

        self.formula: str = None
        self.atoms: List[str] = []
        self.basis_set: str = None
        self.xc_functional: str = None
        self.system: str = None
        self.crystal_system: str = None
        self.spacegroup: str = None
        self.code_name: str = None
        self.code_version: str = None

        for key, value in kwargs.items():
            setattr(self, key, value)

    def to_dict(self):
        return {
            key: value for key, value in self.__dict__.items()
            if value is not None
        }

    def apply_user_metadata(self, metadata: dict):
        """
        Applies a user provided metadata dict to this calc.
        """
        self.pid = metadata.get('_pid')
        self.comment = metadata.get('comment')
        self.upload_time = metadata.get('_upload_time')
        uploader_id = metadata.get('_uploader')
        if uploader_id is not None:
            self.uploader = utils.POPO(id=uploader_id)
        self.references = [utils.POPO(value=ref) for ref in metadata.get('references', [])]
        self.with_embargo = metadata.get('with_embargo', False)
        self.coauthors = [
            utils.POPO(id=user) for user in metadata.get('coauthors', [])]
        self.shared_with = [
            utils.POPO(id=user) for user in metadata.get('shared_with', [])]
        self.datasets = [
            utils.POPO(id=ds['id'], doi=utils.POPO(value=ds.get('_doi')), name=ds.get('_name'))
            for ds in metadata.get('datasets', [])]
