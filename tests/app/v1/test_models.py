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
from typing import List
import pytest
from pydantic import BaseModel, Field, ValidationError
import yaml
import sys

from nomad.utils import strip
from nomad.app.v1.models.graph import GraphRequest
from nomad.app.v1.models.graph.utils import (
    generate_request_model,
    mapped
)


@pytest.fixture()
def path_ref_prefix(monkeypatch):
    monkeypatch.setattr("nomad.app.v1.models.graph.utils.ref_prefix", "#/definitions")


def assert_path(data: BaseModel, path: str):
    if ':' in path:
        path, type_or_value = path.split(':', 1)
    failed_message = f"Cannot find {path}."
    segments = path.split(".")
    current = data
    for segment in segments:
        if isinstance(current, BaseModel):
            value = getattr(current, segment, None)
        elif isinstance(current, dict):
            value = current.get(segment, None)
        else:
            value = current

        assert value is not None, f"{failed_message} Key {segment} in {current} does not exist."

        current = value

    if isinstance(current, BaseModel):
        type_name = current.__class__.__name__
        assert type_name == type_or_value, f'{failed_message} Wrong type; got {type_name} expected {type_or_value}'
    else:
        assert str(current) == type_or_value, f'{failed_message} Wrong value; got {str(current)} expected {type_or_value}'


def test_module():
    for model_name in ['GraphRequest', 'UploadsRequest', 'UploadRequest', 'UploadsResponse']:
        assert hasattr(sys.modules['nomad.app.v1.models.graph.graph_models'], model_name)


@pytest.mark.parametrize(
    "request_yaml, paths, error_path",
    [
        pytest.param(
            """
            users:
                me:
                    username: '*'
                    uploads:
                        m_request:
                            query:
                                upload_name: ['test']
                                is_processing: true
                                is_published: false
                            pagination:
                                page_size: 23
                                page: 1
                        '*':
                            upload_name: '*'
                            files:
                                foo:
                                    bar:
                                        m_is: File
                                        path: '*'
                                        size: '*'
                                        entry:
                                            entry_id: '*'
                            entries:
                                '*':
                                    archive:
                                        data:
                                            my_ref:
                                                m_is: MValue
                                                ref_value:
                                                    '*': '*'
                    datasets:
                        m_request:
                            query:
                                dataset_name: 'My dataset'
                        '*':
                            doi: '*'
            search:
                m_request:
                    owner: visible
                    query:
                        results.material.elements: Fe
                '*':
                    entry:
                        process_status: '*'
                    '*': '*'
            """,
            [
                "users.m_children.me.uploads.m_request.query.is_published:False",
                "users.m_children.me.uploads.m_request.pagination.page_size:23",
                "users.m_children.me.uploads.m_children.*.upload_name:*",
                "users.m_children.me.uploads.m_children.*.entries.m_children.*.archive.m_children.data.m_children.my_ref:MValueRequest",
                "users.m_children.me.uploads.m_children.*.entries.m_children.*.archive.m_children.data.m_children.my_ref.ref_value.m_children.*:*",
                "users.m_children.me.uploads.m_children.*.files.m_children.foo.m_children.bar.entry.entry_id:*",
                "users.m_children.me.datasets.m_request.query:DatasetQuery",
                "users.m_children.me.datasets.m_children.*.doi:*",
                "search.m_children.*.entry.process_status:*",
                "search.m_children.*.entry:EntryRequest",
                "search.m_children.*.*:*"
            ],
            None,
            id="ok",
        ),
        pytest.param(
            """
                does_not_exist:
            """,
            [],
            "does_not_exist",
            id="extra-not-allowed",
        ),
        pytest.param(
            """
                uploads:
                    '*':
                        upload_name: 'wrong type'
            """,
            [],
            "uploads.*.upload_name",
            id="only-*-allowed-str",
        ),
        pytest.param(
            """
                uploads:
                    '*':
                        n_entries: 'wrong type'
            """,
            [],
            "uploads.*.n_entries",
            id="only-*-allowed-other",
        ),
    ],
)
def test_validation(request_yaml: str, paths: List[str], error_path: str):
    try:
        request = GraphRequest.parse_obj(yaml.safe_load(strip(request_yaml)))
    except ValidationError as error:
        assert error_path, str(error)
        assert len(error.errors()) == 1
        loc = [loc for loc in error.errors()[0]["loc"] if loc != "__root__"]
        assert loc == error_path.split(".")
    else:
        assert (
            not error_path
        ), f"Expected validation error in {error_path}, but data passed validation."
        export_kwargs = dict(
            exclude_unset=True, exclude_defaults=False, exclude_none=False
        )
        print(request.json(indent=2, **export_kwargs))
        for path in paths:
            assert_path(request, path)


def test_mapped(path_ref_prefix):
    class MyBase(BaseModel):
        p1: str

    class MySource(MyBase):
        p2: str = Field("p2", description="docs")
        p4: int = 0

    class MyTarget(mapped(MySource, p2="p3", p4=str)):
        p2: int

    target = MyTarget(p1="1", p2=2, p4="p4")
    assert target.p1 == "1"
    assert target.p2 == 2
    assert target.p3 == "p2"
    assert target.p4 == "p4"
    assert MyTarget.__fields__["p3"].field_info.description == "docs"


class Recursive(BaseModel):
    m_children: Recursive


def test_recursive_model(path_ref_prefix):
    root_model = generate_request_model(Recursive)

    root_schema = root_model.schema()
    path_request_schema = root_schema["definitions"]["RecursiveRequest"]
    assert (
        path_request_schema["additionalProperties"]["anyOf"][0]["$ref"]
        == "#/definitions/RecursiveRequest"
    )

    root_model.parse_obj({"m_children": {"name": {"m_children": {}}}})


def test_request_model(path_ref_prefix):
    from nomad.app.v1.models.graph import GraphRequest as root_model

    assert root_model.__module__ == "nomad.app.v1.models.graph.graph_models"
    assert root_model.__name__ == "GraphRequest"

    root_schema = root_model.schema()
    defs = root_schema["definitions"]
    assert "UploadRequest" in defs
    assert [
        type["$ref"]
        for type in defs["UploadsRequest"]["additionalProperties"]["anyOf"]
        if "$ref" in type
    ] == [
        "#/definitions/UploadRequest",
        "#/definitions/UploadRequestOptions",
    ]
    assert "required" not in defs["UploadsRequest"]

    root_model.parse_obj(
        yaml.safe_load(
            strip(
                """
                    uploads:
                        m_request:
                            pagination:
                                page: 1
                        '*':
                            upload_id: '*'
                """
            )
        )
    )


def test_response_model(path_ref_prefix):
    from nomad.app.v1.models.graph.graph_models import GraphResponse as root_model

    assert root_model.__module__ == "nomad.app.v1.models.graph.graph_models"
    assert root_model.__name__ == "GraphResponse"

    root_schema = root_model.schema()
    defs = root_schema["definitions"]
    assert "UploadResponse" in defs
    assert [
        type["$ref"]
        for type in defs["UploadsResponse"]["additionalProperties"]["anyOf"]
        if "$ref" in type
    ] == [
        f"#/definitions/UploadResponse",
        f"#/definitions/UploadResponseOptions",
    ]
    # assert defs["UploadsResponse"]["required"] == ["m_response"]

    root_model.parse_obj(
        yaml.safe_load(
            strip(
                """
                    uploads:
                        m_response:
                            query: {}
                            pagination:
                                total: 1
                                page: 1
                                page_offset: 0
                        '1':
                            upload_id: '1'
                """
            )
        )
    )
