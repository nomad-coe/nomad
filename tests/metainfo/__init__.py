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

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class MTypes:
    # todo: account for bytes which cannot be naturally serialized to JSON
    primitive = {
        str: lambda v: None if v is None else str(v),
        int: lambda v: None if v is None else int(v),
        float: lambda v: None if v is None else float(v),
        complex: lambda v: None if v is None else complex(v),
        bool: lambda v: None if v is None else bool(v),
        np.bool_: lambda v: None if v is None else bool(v),
    }

    primitive_name = {v.__name__: v for v in primitive} | {
        'string': str,
        'boolean': bool,
    }

    int_numpy = {
        np.int8,
        np.int16,
        np.int32,
        np.int64,
        np.uint8,
        np.uint16,
        np.uint32,
        np.uint64,
    }
    int_python = {int}
    int = int_python | int_numpy
    float_numpy = {np.float16, np.float32, np.float64}
    complex_numpy = {np.complex64, np.complex128}
    float_python = {float}
    complex_python = {complex}
    float = float_python | float_numpy
    complex = complex_python | complex_numpy
    num_numpy = int_numpy | float_numpy | complex_numpy
    num_python = int_python | float_python | complex_python
    num = num_python | num_numpy
    str_numpy = {np.str_}
    bool_numpy = {np.bool_}
    bool = {bool, np.bool_}
    numpy = num_numpy | str_numpy | bool_numpy
    str = {str} | str_numpy
