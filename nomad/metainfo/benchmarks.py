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

from nomad.metainfo import MSection, Quantity


class Section(MSection):
    foo = Quantity(name='foo', type=int)


bar_section = Section()
bar_section.foo = 1


class Regular:
    def __init__(self):
        self.foo = 1


obj = Regular()


class Property:
    def __init__(self):
        self.foo = 1

    @property
    def foo(self,):
        return self._foo

    @foo.setter
    def foo(self, value):
        self._foo = value


prop = Property()

dct = dict(foo=1)


def mi(iterations=10):
    bar_section.foo = 1
    for _ in range(0, iterations):
        bar_section.foo += bar_section.foo


def python_obj(iterations=10):
    obj.foo = 1
    for _ in range(0, iterations):
        obj.foo += obj.foo


def python_dct(iterations=10):
    dct['foo'] = 1
    for _ in range(0, iterations):
        dct['foo'] += dct['foo']


def python_property(iterations=10):
    prop.foo = 1
    for _ in range(0, iterations):
        prop.foo += prop.foo


def test_mi(benchmark):
    benchmark(mi)


def test_python_obj(benchmark):
    benchmark(python_obj)


def test_python_dct(benchmark):
    benchmark(python_dct)


def test_python_prop(benchmark):
    benchmark(python_property)
