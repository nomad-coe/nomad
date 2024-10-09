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
from typing import Optional, List, Union, Any, Literal
from pydantic import BaseModel, Field, Extra

from nomad.metainfo import Package
from nomad.graph.model import (
    RequestConfig,
    DatasetQuery,
    MetainfoQuery,
    MetainfoPagination,
)
from nomad.metainfo.pydantic_extension import PydanticModel
from nomad.datamodel.data import User as UserModel
from nomad.app.v1.models.models import Metadata, MetadataResponse
from nomad.app.v1.routers.datasets import Dataset as DatasetV1, DatasetPagination
from nomad.app.v1.routers.uploads import (
    UploadProcData,
    UploadProcDataPagination,
    UploadProcDataQuery,
    PaginationResponse,
    EntryProcData,
    EntryProcDataPagination,
)

from nomad.app.v1.models.graph.utils import (
    generate_request_model,
    generate_response_model,
    mapped,
)


class Error(BaseModel):
    error_type: str
    message: str


RecursionOptions = RequestConfig

DirectoryRequestOptions = RequestConfig


class DirectoryResponseOptions(BaseModel):
    pagination: PaginationResponse


class GraphDirectory(BaseModel):
    m_errors: List[Error]
    m_is: Literal['Directory']
    m_request: DirectoryRequestOptions
    m_response: DirectoryResponseOptions
    m_children: Union[GraphDirectory, GraphFile]


class GraphFile(BaseModel):
    m_errors: List[Error]
    m_is: Literal['File']
    m_request: DirectoryRequestOptions
    path: str
    size: int
    entry: Optional[GraphEntry]
    # The old API also had those, but they can be grabbed from entry:
    # parser_name, entry_id, archive
    # This is similar to the question for "m_parent" in Directory. At least we need
    # to navigate from Entry to mainfile to directory, but we could also but a
    # mainfile_directory into Entry?
    parent: GraphDirectory


class MSection(BaseModel):
    m_errors: List[Error]
    m_request: RecursionOptions
    m_def: MDef
    m_children: Any


class MDef(MSection):
    m_def: str  # type: ignore
    m_def_id: str


class GraphEntry(mapped(EntryProcData, mainfile='mainfile_path', entry_metadata=None)):  # type: ignore
    m_errors: List[Error]
    mainfile: GraphFile
    upload: GraphUpload
    archive: MSection
    metadata: GraphEntryMetadata


class EntriesRequestOptions(BaseModel):
    # The old API does not support any queries
    pagination: Optional[EntryProcDataPagination]


class EntriesResponseOptions(BaseModel):
    pagination: Optional[PaginationResponse]
    # The "upload" was only necessary, because in the old API you would not get the upload.
    # In the graph API, the upload would be the parent anyways
    # upload: Upload


class GraphEntries(BaseModel):
    m_request: EntriesRequestOptions
    m_response: EntriesResponseOptions
    m_children: GraphEntry


class GraphUser(
    UserModel.m_def.m_get_annotation(PydanticModel).model,  # type: ignore
    extra=Extra.forbid,
):
    # This is more complicated as the user can have different roles in different uploads.
    # This would only refer to uploads with the user as main_author.
    # For many clients and use-cases uploads.m_request.query will be the
    # more generic or only option
    uploads: Optional[GraphUploads]
    datasets: Optional[GraphDatasets]


class GraphUsers(BaseModel):
    m_errors: List[Error]
    m_children: GraphUser


class GraphUpload(
    mapped(  # type: ignore
        UploadProcData,
        entries='n_entries',
        main_author=GraphUser,
        coauthors=List[GraphUser],
        reviewers=List[GraphUser],
        viewers=List[GraphUser],
        writers=List[GraphUser],
    ),
    extra=Extra.forbid,
):
    # The old API includes some extra data here:
    processing_successful: int = Field(
        description='Number of entries that has been processed successfully.'
    )
    processing_failed: int = Field(
        description='Number of entries that failed to process.'
    )

    entries: GraphEntries = Field(description='The entries contained in this upload.')
    files: GraphDirectory = Field(
        description="This upload's root directory for all files (raw data)."
    )


class UploadRequestOptions(BaseModel):
    pagination: Optional[UploadProcDataPagination]
    query: Optional[UploadProcDataQuery]


class UploadResponseOptions(BaseModel):
    pagination: Optional[PaginationResponse]
    query: Optional[UploadProcDataQuery]


class GraphUploads(BaseModel):
    m_request: UploadRequestOptions
    m_response: UploadResponseOptions
    m_errors: List[Error]
    m_children: GraphUpload


class GraphEntryMetadata(BaseModel, extra=Extra.allow):
    entry: GraphEntry


class SearchRequestOptions(BaseModel):
    query: Optional[Metadata]


class SearchResponseOptions(BaseModel):
    query: Optional[MetadataResponse]


class GraphSearch(BaseModel):
    m_request: SearchRequestOptions
    m_response: SearchResponseOptions
    m_errors: List[Error]
    m_children: GraphEntryMetadata


class GraphDataset(mapped(DatasetV1, query=None, entries=None)):  # type: ignore
    pass


class DatasetRequestOptions(BaseModel):
    pagination: Optional[DatasetPagination]
    query: Optional[DatasetQuery]


class DatasetResponseOptions(BaseModel):
    pagination: Optional[PaginationResponse]
    query: Optional[DatasetQuery]


class GraphDatasets(BaseModel):
    m_request: DatasetRequestOptions
    m_response: DatasetResponseOptions
    m_errors: List[Error]
    m_children: GraphDataset


class MetainfoRequestOptions(BaseModel):
    pagination: Optional[MetainfoPagination]
    query: Optional[MetainfoQuery]


class MetainfoResponseOptions(BaseModel):
    pagination: Optional[PaginationResponse]
    query: Optional[MetainfoQuery]


class GraphMetainfo(BaseModel):
    m_request: MetainfoRequestOptions
    m_response: MetainfoResponseOptions
    m_errors: List[Error]
    m_children: MSection


class Graph(BaseModel):
    users: GraphUsers
    entries: GraphEntries
    uploads: GraphUploads
    datasets: GraphDatasets
    search: GraphSearch
    metainfo: GraphMetainfo


GraphRequest = generate_request_model(Graph)
GraphResponse = generate_response_model(Graph)
