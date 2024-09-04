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
