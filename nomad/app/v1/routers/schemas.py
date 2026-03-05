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
import importlib
from enum import Enum
from typing import Annotated

from fastapi import APIRouter, HTTPException, Path, Query, status
from fastapi.responses import JSONResponse

router = APIRouter()


class APITag(str, Enum):
    DEFAULT = 'schemas'


class SerializationFormat(str, Enum):
    JSONSCHEMA = 'jsonschema'
    METAINFO = 'metainfo'


@router.get(
    '/{schema_id}',
    tags=[APITag.DEFAULT],
    summary='Return a serialization of a specific data schema.',
)
async def get_schema(
    schema_id: Annotated[
        str,
        Path(
            description="""
Schema identifier. For now, we only support using the qualified name of a section,
e.g. `package_name.schema_packages.calculations.MySchema`. This will return the latest
schema registered into the system.
"""
        ),
    ],
    format: Annotated[
        SerializationFormat,
        Query(
            description="""
Format for the returned schema. Available formats:
  - `jsonschema` (default) - JSON Schema (Draft 2020-12) representation.

    Always-included keywords:
      - `"$schema"`: The JSONSchema version
      - `"$id"`: Identifier as a resolvable URL

    Optional extras (added only when present):
      - `"title"`: Name of the quantity/section
      - `"description"`: Description of the quantity/section
      - `"unit"`: Unit of a quantity as string

    Type Data Keywords:
      - `"type"`                   - type of the quantity/section.
      - Scalar quantities are mapped to the corresponding JSON type, e.g. `number`.
      - Arrays: `shape` is walked from left to right, wrapping each dimension
        in a nested ``{"type": "array", …, "items": {…}}``. Array shape is
        mapped into `minItems`/`maxItems` as follows:
          - `n` (int)              - `minItems = maxItems = n`
          - `"*"`                  - unbounded: no min/max keys
          - `"a.."`                - `minItems = a`
          - `"..b"`                - `maxItems = b`
          - `"a..b"`               - `minItems = a`, `maxItems = b`
          - any other string       - treated as `"*"`
      - Optional fields for quantities:
          - `"default"`      - Default value of a quantity
          - `"enum"`         - List of allowed values for an Enum quantity
          - `"minimum"`      - Minimum value for a quantity (from ELN annotation)
          - `"maximum"`      - Maximum value for a quantity (from ELN annotation)
      - Reference quantities are mapped to type strings
      - Sections have the `object` type
  - `metainfo` - The NOMAD-specific serialization format with full fidelity to the original schema.
"""
        ),
    ] = SerializationFormat.JSONSCHEMA,
):
    """
    Returns the serialized scehma for the given id. The returned
    schema is serialized in the format specified by the `format` query parameter.
    """

    def resolve_m_def(m_def: str):
        """
        Resolve Section from qualified name (m_def) such as:
            package_name.schema_packages.calculations.MySchema
        """
        parts: list[str] = m_def.split('.')
        module_path: str = '.'.join(parts[:-1])
        class_name: str = parts[-1]

        try:
            module = importlib.import_module(module_path)
            section = getattr(module, class_name)
        except (ImportError, AttributeError, ValueError) as e:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f'Class {class_name} does not exist in module {module_path}.',
            ) from e

        if not hasattr(section, 'm_def'):
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f'{section=} does not have metainfo definition.',
            )

        return section

    section = resolve_m_def(m_def=schema_id)

    if format == SerializationFormat.JSONSCHEMA:
        return JSONResponse(
            content=section.m_def.m_to_json_schema(),
            media_type='application/schema+json',
        )
    elif format == SerializationFormat.METAINFO:
        return JSONResponse(
            content=section.m_def.m_to_dict(),
            media_type='application/json',
        )
