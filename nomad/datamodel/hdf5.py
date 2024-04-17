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
from typing import Any, cast
import h5py
import re

from nomad.utils import get_logger
from nomad.metainfo import DataType, MSection, Quantity
from nomad.datamodel.context import ServerContext, Context
from nomad.datamodel import EntryArchive

LOGGER = get_logger(__name__)


def match_hdf5_reference(reference: str):
    """
    Match reference to HDF5 upload path syntax.
    """
    return re.match(
        r'(?:.*?/uploads/(?P<upload_id>.+?)/(?P<directory>.+?)/)*(?P<file_id>.+?)#(?P<path>.+)',
        reference,
    )


def read_hdf5_dataset(hdf5_file: h5py.File, path: str) -> h5py.Dataset:
    """
    Read an HDF5 dataset given by path.
    """
    match = match_hdf5_reference(path)
    return (
        hdf5_file[match['file_id']] if match['file_id'] in hdf5_file else hdf5_file
    )[match['path']]


def write_hdf5_dataset(value: Any, hdf5_file: h5py.File, path: str) -> str:
    """
    Write data to HDF5 file.
    """
    segments = path.rsplit('/', 1)
    group = hdf5_file.require_group(segments[0]) if len(segments) == 2 else hdf5_file
    dataset = group.require_dataset(
        segments[-1],
        shape=value.shape if hasattr(value, 'shape') else (),
        dtype=value.dtype if hasattr(value, 'dtype') else None,
    )
    dataset[...] = value.magnitude if hasattr(value, 'magnitude') else value
    return dataset


class _HDF5Reference(DataType):
    @staticmethod
    def _get_upload_files(archive: EntryArchive, path: str):
        match = match_hdf5_reference(path)
        file_id = match['file_id']

        if not file_id:
            LOGGER.error('Invalid HDF5 path.')

        from nomad import files

        upload_id, _ = ServerContext._get_ids(archive, required=True)
        return match, files.UploadFiles.get(upload_id)

    @staticmethod
    def write_dataset(archive: EntryArchive, value: Any, path: str) -> None:
        """
        Write value to HDF5 file specified in path following the form
        filename.h5#/path/to/dataset. upload_id is resolved from archive.
        """

        match, upload_files = HDF5Reference._get_upload_files(archive, path)
        mode = 'r+b' if upload_files.raw_path_is_file(match['file_id']) else 'wb'
        with h5py.File(upload_files.raw_file(match['file_id'], mode), 'a') as f:
            write_hdf5_dataset(value, f, match['path'])

    @staticmethod
    def read_dataset(archive: EntryArchive, path: str) -> Any:
        """
        Read HDF5 dataset from file specified in path following the form
        filename.h5#/path/to/dataset. upload_id is resolved from archive.
        """
        match, upload_files = HDF5Reference._get_upload_files(archive, path)
        with h5py.File(upload_files.raw_file(match['file_id'], 'rb')) as f:
            return read_hdf5_dataset(f, path)[()]


class _HDF5Dataset(DataType):
    def serialize(self, section: MSection, quantity_def: Quantity, value: Any) -> str:
        if isinstance(value, h5py.Dataset):
            return f'{value.file.filename}#{value.name}'

        section_root: MSection = section.m_root()
        if not section_root.m_context:
            LOGGER.error(
                'Cannot serialize HDF5 value.',
                data=dict(quantity_definition=quantity_def),
            )
            return None

        context = cast(Context, section_root.m_context)

        path = f'{section.m_path()}/{quantity_def.name}'
        upload_id, entry_id = ServerContext._get_ids(section.m_root(), required=True)

        with context.open_hdf5_file(section, quantity_def, value, 'w') as file_object:
            hdf5_file = h5py.File(file_object, 'a')
            try:
                write_hdf5_dataset(value, hdf5_file, path)
                return f'/uploads/{upload_id}/archive/{entry_id}#{path}'
            except Exception as e:
                LOGGER.error('Cannot write HDF5 dataset', exc_info=e)
                return None
            finally:
                hdf5_file.close()

    def deserialize(
        self, section: MSection, quantity_def: Quantity, value: str
    ) -> h5py.Dataset:
        match = match_hdf5_reference(value)
        if not match:
            return None

        section_root: MSection = section.m_root()
        if not section_root.m_context:
            LOGGER.error(
                'Cannot resolve HDF5 reference.',
                data=dict(quantity_definition=quantity_def),
            )
            return None

        context = cast(Context, section_root.m_context)

        file_object = context.open_hdf5_file(section, quantity_def, value, 'r')
        hdf5_file = context.file_handles.get(file_object.name)
        if not hdf5_file:
            hdf5_file = h5py.File(file_object)
            context.file_handles[file_object.name] = hdf5_file

        return read_hdf5_dataset(hdf5_file, value)


HDF5Reference = _HDF5Reference()
HDF5Dataset = _HDF5Dataset()
