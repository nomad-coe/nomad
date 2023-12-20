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
from typing import (
    Dict,
    List,
    Optional,
    Type,
    Literal,
    Union,
    Any,
    Callable,
    ForwardRef,
    get_type_hints,
    get_origin,
    get_args,
    cast,
)
from datetime import datetime
from pydantic import (
    BaseModel,
    BaseConfig,
    create_model,
    Extra,
    Field,
    root_validator,
    ValidationError,
    parse_obj_as,
)
from pydantic.error_wrappers import ErrorWrapper
from pydantic.typing import evaluate_forwardref
from pydantic.config import inherit_config
import sys


ref_prefix = "#/components/schemas"
request_suffix = "Request"
response_suffix = "Response"


class _DictModel(BaseModel):
    @classmethod
    def process_extra(cls, values):
        m_children = values.setdefault("m_children", {})
        type_ = cls.__fields__["m_children"].type_
        for name in list(values):
            if name not in cls.__fields__:
                value = values[name]
                values.pop(name)
                try:
                    m_children[name] = parse_obj_as(type_, value)
                except ValidationError as exc:
                    # m_children is always a Union and the last possible type is
                    # Literal['*']. Respectively the last validation errors comes from
                    # this type. It is usually confusing and not helpful to the user.
                    # Therefore, we pop it.
                    if len(exc.raw_errors) > 1:
                        exc.raw_errors.pop()  # pylint: disable=no-member
                    raise ValidationError([ErrorWrapper(exc, loc=name)], cls)

        return values

    class Config:
        extra = Extra.allow

        @staticmethod
        def schema_extra(schema: dict[str, Any], model: Type[_DictModel]) -> None:
            if "m_children" not in model.__annotations__:
                raise TypeError(
                    f"No m_children field defined for dict model {model.__name__}. "
                )
            children_annotation = model.__annotations__["m_children"]
            value_type = get_args(get_args(children_annotation)[0])[1]
            if value_type is None:
                raise TypeError(
                    f"Could not determine m_children's type. Did you miss to call update_forward_refs()?"
                )

            if get_origin(value_type) == Union:
                value_types = get_args(value_type)
            else:
                value_types = (value_type,)

            types = []
            for value_type in value_types:
                if isinstance(value_type, ForwardRef):
                    value_type = value_type.__forward_value__

                if value_type == Literal["*"]:
                    types.append({"enum": ["*"], "type": "string"})
                else:
                    types.append({"$ref": f"{ref_prefix}/{value_type.__name__}"})

            if "properties" in schema:
                for property in schema["properties"].values():
                    if "$ref" in property:
                        types.append(property)

                schema["properties"].pop("m_children")

            schema["additionalProperties"] = {"anyOf": types}


def _get_request_type(type_hint: Any, ns: ModelNamespace) -> Any:
    origin, args = get_origin(type_hint), get_args(type_hint)

    if origin is list or type_hint in [str, int, bool, datetime, Any]:
        return Literal["*"]

    if origin is None and issubclass(type_hint, BaseModel):
        return _generate_model(type_hint, request_suffix, _get_request_type, ns)

    if origin is dict:
        key_type, value_type = args
        return Dict[key_type, _get_request_type(value_type, ns)]  # type: ignore

    # This is about Optional[T], which is translated to Union[None, T]
    if origin is Union and len(args) == 2 and isinstance(None, args[1]):
        return _get_request_type(args[0], ns)

    if origin is Union:
        union_types = tuple(_get_request_type(type_, ns) for type_ in args)
        return Union[union_types]  # type: ignore

    raise NotImplementedError(type_hint)


def _get_response_type(type_hint: Any, ns: ModelNamespace) -> Any:
    origin, args = get_origin(type_hint), get_args(type_hint)
    if type_hint in [str, int, bool, datetime, Any]:
        return type_hint

    if origin is None and issubclass(type_hint, BaseModel):
        return _generate_model(type_hint, response_suffix, _get_response_type, ns)

    if origin is list:
        value_type = args[0]
        return List[_get_response_type(value_type, ns)]  # type: ignore

    if origin is dict:
        key_type, value_type = args
        # TODO is this really necessary?
        if value_type == type_hint:
            # We have detected direct type recursion, like in
            # Path = Dict[str, 'Path']
            return type_hint
        return Dict[key_type, _get_response_type(value_type, ns)]  # type: ignore

    # This is about Optional[T], which is translated to Union[None, T]
    if origin is Union and len(args) == 2 and isinstance(None, args[1]):
        return _get_response_type(args[0], ns)

    if origin is Union:
        union_types = tuple(_get_response_type(type_, ns) for type_ in args)
        return Union[union_types]  # type: ignore

    raise NotImplementedError(type_hint)


ModelNamespace = Dict[str, Union[Type[BaseModel], ForwardRef]]


def _generate_model(
    source_model: Type[BaseModel],
    suffix: str,
    generate_type: Callable[[type, ModelNamespace], type],
    ns: ModelNamespace,
    **kwargs,
):
    # We need to populate a forward ref for the model in the ns use it in recursion cases.
    result_model_name = f"{source_model.__name__}{suffix}"
    is_ns_origin = len(ns) == 0
    if result_model_name not in ns:
        ns[result_model_name] = ForwardRef(result_model_name)
    else:
        return ns[result_model_name]

    type_hints = get_type_hints(source_model)
    fields = dict(**kwargs)

    for field_name, type_hint in type_hints.items():

        if field_name.startswith("__"):
            continue

        if field_name == "m_children":
            origin, args = get_origin(type_hint), get_args(type_hint)
            if origin is Union:
                types = args
            else:
                types = (type_hint,)
            if not all(isinstance(type_, type) and issubclass(type_, BaseModel) for type_ in types):
                raise TypeError(
                    "Only Pydantic model classes (or Unions thereof) are supported as m_children types."
                )
            value_types = tuple(_generate_model(type_, suffix, generate_type, ns) for type_ in types)
            # TODO we always add Literal['*'] at the end. Maybe it should be configurable
            # which models want to support '*' values for their children?
            value_type = Union[value_types + (Literal['*'],)]  # type: ignore
            fields["m_children"] = (Optional[Dict[str, cast(Type, value_type)]], None)  # type: ignore
            continue

        if field_name == "m_request":
            if suffix == request_suffix:
                fields[field_name] = (Optional[type_hint], None)
            continue

        if field_name == "m_response":
            if suffix == response_suffix:
                fields[field_name] = (Optional[type_hint], None)
            continue

        if field_name == 'm_is':
            fields[field_name] = (Optional[type_hint], None)
            continue

        if field_name == 'm_errors':
            if suffix == response_suffix:
                fields[field_name] = (Optional[Union[type_hint]], None)  # type: ignore
            continue

        if field_name.startswith('m_') and field_name not in ['m_def']:
            raise NotImplementedError(f'The internal field {field_name} is not implemented.')

        fields[field_name] = (Optional[generate_type(type_hint, ns)], None)

    config = source_model.__config__
    if config.extra == Extra.ignore and 'm_children' not in fields:
        config = inherit_config(
            type('Config', (BaseConfig,), dict(extra=Extra.forbid)), config
        )

    validators = {}
    if 'm_children' in fields:
        config = inherit_config(_DictModel.__config__, config)
        if suffix == request_suffix:
            validators = {
                'process_extra': root_validator(  # type: ignore
                    _DictModel.process_extra.__func__, # type: ignore
                    pre=True,
                    allow_reuse=True,
                )
            }

    result_model = create_model(
        result_model_name,
        __module__=source_model.__module__,
        __validators__=validators,
        __config__=config,
        **fields,
    )

    # We need to replace the forward ref in the ns with the finished model. We also
    # need to update all forward refs after the whole model has been created.
    ns[result_model_name] = result_model
    if is_ns_origin:
        for model in ns.values():
            if isinstance(model, type):
                model.update_forward_refs(**ns)
                # There is a bug in pydantics BaseModel.update_forward_refs and it does not
                # recognize forward refs in Union types. Therefore we do our own impl.
                # https://github.com/pydantic/pydantic/issues/3345
                for field in model.__fields__.values():
                    if get_origin(field.type_) is Union:
                        union_types = tuple(
                            evaluate_forwardref(type_, {}, ns) if type_.__class__ == ForwardRef else type_
                            for type_ in get_args(field.type_)
                        )
                        field.type_ = Union[union_types]  # type: ignore

    assert getattr(sys.modules[source_model.__module__], result_model_name, result_model) == result_model, \
        f'Model class with name {result_model_name} already exists.'
    setattr(sys.modules[source_model.__module__], result_model_name, result_model)

    return result_model


def mapped(
        model: Type[BaseModel], **mapping: Union[str, type]
) -> Type[BaseModel]:
    """
    Creates a new pydantic model based on the given model. The mapping argument allows
    to either change the name of a field in the input model or change the type of a field
    in the given input model or remove the field.

    Arguments:
        model: the input model.
        **kwargs: field names and either the new field name or the new field type or
            None to remove the field.

    Returns:
        a pydantic model with the mapped field and the same base as the input model
    """

    def create_field(field_info):
        return Field(
            default=field_info.default,
            alias=field_info.alias,
            title=field_info.title,
            description=field_info.description,
        )

    fields = {}
    for name, field in model.__fields__.items():
        if name not in mapping:
            fields[name] = (field.type_, create_field(field.field_info))
            continue

        new_name_or_type_ = mapping[name]
        old_field = model.__fields__[name]

        if new_name_or_type_ is None:
            continue

        if isinstance(new_name_or_type_, str):
            new_name = new_name_or_type_
            type_ = old_field.type_
        else:
            new_name = name
            type_ = new_name_or_type_

        fields[new_name] = (type_, create_field(old_field.field_info))

    return create_model(  # type: ignore
        model.__name__, **fields, __module__=model.__module__, __base__=model.__base__
    )


def generate_request_model(source_model: Type[BaseModel]):
    return _generate_model(source_model, request_suffix, _get_request_type, dict())


def generate_response_model(source_model: Type[BaseModel]):
    return _generate_model(source_model, response_suffix, _get_response_type, dict())
