"""Make sure file has the NOMAD copyright header."""

import argparse
from pathlib import Path

HEADER_BODY = """#
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
"""

EXPECTED_HEADER = HEADER_BODY + '\n'


def _split_prefix(lines: list[str]) -> tuple[list[str], list[str]]:
    prefix = []

    if lines and lines[0].startswith('#!'):
        prefix.append(lines.pop(0))
    if lines and 'coding' in lines[0]:
        prefix.append(lines.pop(0))

    return prefix, lines


def _normalize_header(path: Path) -> None:
    text = path.read_text(encoding='utf-8')
    if text.startswith(EXPECTED_HEADER):
        return

    lines = text.splitlines(keepends=True)
    prefix, lines = _split_prefix(lines)
    body = ''.join(lines)

    if body.startswith(HEADER_BODY):
        remainder = body[len(HEADER_BODY) :].lstrip('\n')
        path.write_text(''.join(prefix) + EXPECTED_HEADER + remainder, encoding='utf-8')
        return

    path.write_text(''.join(prefix) + EXPECTED_HEADER + body, encoding='utf-8')


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--fix',
        action='store_true',
        help='Automatically fix missing or malformed headers in-place.',
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    root = Path(__file__).resolve().parents[1]
    missing_header = []

    for directory in (root / 'nomad', root / 'tests'):
        for path in sorted(directory.rglob('*.py')):
            if args.fix:
                _normalize_header(path)
            if not path.read_text(encoding='utf-8').startswith(EXPECTED_HEADER):
                missing_header.append(path.relative_to(root).as_posix())

    if not missing_header:
        return 0

    print('These Python files are missing the NOMAD copyright header.')
    print('Run locally with --fix to prepend it automatically:')
    print('python scripts/check_copyright_headers.py --fix')
    print()
    print('\n'.join(missing_header))
    return 1


if __name__ == '__main__':
    raise SystemExit(main())
