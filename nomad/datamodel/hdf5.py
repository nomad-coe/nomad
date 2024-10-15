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
from __future__ import annotations

from typing import Any
import h5py
import re

import numpy as np
import pint
from h5py import File

from nomad.metainfo.data_type import NonPrimitive
from nomad.utils import get_logger

LOGGER = get_logger(__name__)


_H5_FILE_ = re.compile(
    r'(?:.*?/uploads/(?P<upload_id>.+?)/(?P<directory>.+?)/)*(?P<file_id>.+?)#(?P<path>.+)'
)


def match_hdf5_reference(reference: str):
    """
    Match reference to HDF5 upload path syntax.
    """
    if not (match := _H5_FILE_.match(reference)):
        return None

    return match.groupdict()


class HDF5Wrapper:
    def __init__(self, file: str, path: str):
        self.file: str = file
        self.path: str = path
        self.handler: h5py.File | None = None

    def __enter__(self):
        self._close()
        self.handler = h5py.File(self.file, 'a')
        return self.handler[self.path]

    def __exit__(self, exc_type, exc_value, traceback):
        self._close()

    def _close(self):
        if self.handler:
            self.handler.close()


class HDF5Reference(NonPrimitive):
    @staticmethod
    def _get_upload_files(archive, path: str):
        match = match_hdf5_reference(path)
        file_id = match['file_id']

        if not file_id:
            LOGGER.error('Invalid HDF5 path.')

        from nomad import files
        from nomad.datamodel.context import ServerContext

        upload_id, _ = ServerContext._get_ids(archive, required=True)
        return file_id, match['path'], files.UploadFiles.get(upload_id)

    @staticmethod
    def write_dataset(archive, value: Any, path: str) -> None:
        """
        Write value to HDF5 file specified in path following the form
        filename.h5#/path/to/dataset. upload_id is resolved from archive.
        """

        file, path, upload_files = HDF5Reference._get_upload_files(archive, path)
        mode = 'r+b' if upload_files.raw_path_is_file(file) else 'wb'
        with h5py.File(upload_files.raw_file(file, mode), 'a') as f:
            f.require_dataset(
                path,
                shape=getattr(value, 'shape', ()),
                dtype=getattr(value, 'dtype', None),
            )[...] = getattr(value, 'magnitude', value)

    @staticmethod
    def read_dataset(archive, path: str) -> Any:
        """
        Read HDF5 dataset from file specified in path following the form
        filename.h5#/path/to/dataset. upload_id is resolved from archive.
        """
        file, path, upload_files = HDF5Reference._get_upload_files(archive, path)
        with h5py.File(upload_files.raw_file(file, 'rb')) as f:
            return (f[file] if file in f else f)[path][()]

    def _normalize_impl(self, value, **kwargs):
        if match_hdf5_reference(value) is not None:
            return value

        raise ValueError(f'Invalid HDF5 reference: {value}.')


class HDF5Dataset(NonPrimitive):
    def _serialize_impl(self, value, **kwargs):
        if isinstance(value, str):
            return value

        section = kwargs.get('section')

        if not (section_context := section.m_root().m_context):
            raise ValueError('Cannot normalize HDF5 value without context.')

        upload_id, entry_id = section_context._get_ids(section.m_root(), required=True)

        return f'/uploads/{upload_id}/archive/{entry_id}#{section.m_path()}/{self._definition.name}'

    def _normalize_impl(self, value, **kwargs):
        """
        In memory, it is represented by either a reference string or a h5py.Dataset.
        The h5py.Dataset handles file access and data storage under the hood for us.
        There is no need to additionally manage file access here.
        """
        section = kwargs.get('section')

        if not (section_context := section.m_root().m_context):
            raise ValueError('Cannot normalize HDF5 value without context.')

        if not isinstance(value, (str, np.ndarray, h5py.Dataset, pint.Quantity)):
            raise ValueError(f'Invalid HDF5 dataset value: {value}.')

        hdf5_path: str = section_context.hdf5_path(section)

        if isinstance(value, str):
            if not (match := match_hdf5_reference(value)):
                # seems to be an illegal reference, do not try to resolve it
                return value

            file, path = match['file_id'], match['path']

            with File(hdf5_path, 'a') as hdf5_file:
                if file in hdf5_file:
                    segment = f'{file}/{path}'
                else:
                    segment = path
        else:
            if isinstance(value, pint.Quantity):
                if self._definition.unit is not None:
                    value = value.to(self._definition.unit).magnitude
                else:
                    value = value.magnitude

            with File(hdf5_path, 'a') as hdf5_file:
                segment = f'{section.m_path()}/{self._definition.name}'
                target_dataset = hdf5_file.require_dataset(
                    segment,
                    shape=getattr(value, 'shape', ()),
                    dtype=getattr(value, 'dtype', None),
                )
                target_dataset[...] = value

        return HDF5Wrapper(hdf5_path, segment)
