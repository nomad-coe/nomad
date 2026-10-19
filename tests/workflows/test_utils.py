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

import pytest

from nomad.workflows.utils import generate_batches


@pytest.mark.parametrize(
    'entry_count, expected_num_batches, expected_max_batch_size',
    [
        (500, 1, 500),
        (5000, 5, 1000),
        (4501, 5, 1000),
        (1_000_000, 1000, 1000),
        (2_000_000, 1000, 2000),
        (1, 1, 1),
        (0, 0, 0),
    ],
)
def test_generate_batches(entry_count, expected_num_batches, expected_max_batch_size):
    entries = list(range(entry_count))
    batches = generate_batches(entries)

    if entry_count == 0:
        assert len(batches) == 0
    else:
        assert len(batches) == expected_num_batches
        assert sum(len(b) for b in batches) == entry_count
        assert max(len(b) for b in batches) <= expected_max_batch_size
