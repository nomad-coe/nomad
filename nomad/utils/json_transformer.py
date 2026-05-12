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
import re
from typing import Any

import jmespath

from nomad.datamodel.metainfo.annotations import Condition, Rule, Rules


class Transformer:
    def __init__(self, mapping_dict: dict[str, 'Rules']):
        self.mapping_dict = mapping_dict

    @staticmethod
    def parse_path(path: str) -> list[str | int]:
        """
        Parses a JMESPath-like path into a list of keys and indices.

        Args:
            path (str): The path string to parse.

        Returns:
            list[Union[str, int]]: A list containing string keys and integer indices.
        """
        pattern = re.compile(r'([^\[\].]+)|\[(\d+)\]')
        parts = []
        for match in pattern.finditer(path):
            key, index = match.groups()
            if key is not None:
                parts.append(key)
            elif index is not None:
                parts.append(int(index))
        return parts

    @staticmethod
    def apply_condition(condition: 'Condition', source: dict[str, Any]) -> bool:
        path = (
            condition.regex_condition.regex_path
            or condition.regex_condition.regex_pattern
        )
        value = jmespath.search(path, source)
        if value is None:
            return False
        value_str = str(value)
        return bool(re.match(condition.regex_condition.regex_pattern, value_str))

    @staticmethod
    def set_value(path: str, value: Any, data: Any) -> None:
        """
        Sets a value in a nested dictionary or list based on the provided path.

        Args:
            path (str): The JMESPath-like path indicating where to set the value.
            value (Any): The value to set.
            data (Any): The target data structure (dict or list).

        Raises:
            TypeError: If there's a mismatch between expected and actual data types.
            IndexError: If a list index is out of bounds.
        """
        parts = Transformer.parse_path(path)
        current = data

        for i, part in enumerate(parts):
            if i == len(parts) - 1:
                if isinstance(part, int):
                    if not isinstance(current, list):
                        raise TypeError(
                            f"Expected list at path '{'.'.join(map(str, parts[:i]))}', "
                            f'but got {type(current).__name__}'
                        )
                    while len(current) <= part:
                        current.append(None)
                    current[part] = value
                else:
                    if not isinstance(current, dict):
                        raise TypeError(
                            f"Expected dict at path '{'.'.join(map(str, parts[:i]))}', "
                            f'but got {type(current).__name__}'
                        )
                    current[part] = value
            else:
                next_part = parts[i + 1]
                if isinstance(part, int):
                    if not isinstance(current, list):
                        raise TypeError(
                            f"Expected list at path '{'.'.join(map(str, parts[:i]))}', "
                            f'but got {type(current).__name__}'
                        )
                    while len(current) <= part:
                        current.append({})
                    if current[part] is None:
                        current[part] = [] if isinstance(next_part, int) else {}
                    current = current[part]
                else:
                    if not isinstance(current, dict):
                        raise TypeError(
                            f"Expected dict at path '{'.'.join(map(str, parts[:i]))}', "
                            f'but got {type(current).__name__}'
                        )
                    if part not in current or current[part] is None:
                        current[part] = [] if isinstance(next_part, int) else {}
                    elif isinstance(next_part, int) and not isinstance(
                        current[part], list
                    ):
                        raise TypeError(
                            f"Expected list at path '{'.'.join(map(str, parts[: i + 1]))}', "
                            f'but got {type(current[part]).__name__}'
                        )
                    elif isinstance(next_part, str) and not isinstance(
                        current[part], dict
                    ):
                        raise TypeError(
                            f"Expected dict at path '{'.'.join(map(str, parts[: i + 1]))}', "
                            f'but got {type(current[part]).__name__}'
                        )
                    current = current[part]

    def resolve_reference(
        self, rule: 'Rule', all_rules: dict[str, 'Rules'], visited=None
    ) -> 'Rule':
        """
        Resolves a rule reference specified in the `use_rule` field.

        Args:
            rule (Rule): The current rule that may reference another rule.
            all_rules (dict[str, Rules]): All available rule sets.
            visited (set, optional): Set of visited references to detect circular dependencies.

        Returns:
            Rule: The resolved rule with overridden fields if a reference exists.
        """
        if visited is None:
            visited = set()

        if rule.use_rule and rule.use_rule.startswith('#'):
            if rule.use_rule in visited:
                raise ValueError(
                    f"Circular reference detected for use_rule '{rule.use_rule}'."
                )
            visited.add(rule.use_rule)

            ref_path = rule.use_rule[1:]
            try:
                mapping_name, rule_name = ref_path.split('.', 1)
            except ValueError:
                raise ValueError(
                    f"Invalid use_rule format: '{rule.use_rule}'. Expected format '#mapping.rule_name'."
                )

            if mapping_name not in all_rules:
                raise ValueError(
                    f"Mapping name '{mapping_name}' not found in the mapping dictionary."
                )

            referenced_rules = all_rules[mapping_name].rules
            if rule_name not in referenced_rules:
                raise ValueError(
                    f"Rule name '{rule_name}' not found in mapping '{mapping_name}'."
                )

            referenced_rule = referenced_rules[rule_name]
            rule = rule.override_fields(referenced_rule)
            rule.use_rule = referenced_rule.use_rule

            rule = self.resolve_reference(rule, all_rules, visited)

        return rule

    def transform_dict(
        self,
        rule: 'Rule',
        source: dict[str, Any],
        target: Any,
        all_rules: dict[str, 'Rules'],
        parent_source_path: str = '',
        parent_target_path: str = '',
        visited=None,
    ) -> Any:
        """
        Transforms the source dictionary into the target based on the provided rule.

        Args:
            rule (Rule): The transformation rule to apply.
            source (dict[str, Any]): The source data.
            target (Any): The target data structure.
            all_rules (dict[str, 'Rules']): All available rule sets.
            parent_source_path (str, optional): The parent source path.
            parent_target_path (str, optional): The parent target path.
            visited (set, optional): Set of visited references to detect circular dependencies.

        Returns:
            Any: The updated target data structure.
        """
        resolved_rule = self.resolve_reference(rule, all_rules, visited)

        current_source_path = resolved_rule.source or parent_source_path
        current_target_path = resolved_rule.target or parent_target_path

        conditions_met = True
        if resolved_rule.conditions:
            conditions_met = all(
                self.apply_condition(cond, source) for cond in resolved_rule.conditions
            )

        source_value = None
        if current_source_path:
            source_value = jmespath.search(current_source_path, source)

        if conditions_met:
            if source_value is not None:
                self.set_value(current_target_path, source_value, target)
            elif resolved_rule.default_value is not None:
                self.set_value(current_target_path, resolved_rule.default_value, target)
        else:
            if resolved_rule.default_value is not None:
                self.set_value(current_target_path, resolved_rule.default_value, target)

        return target

    def dict_to_dict(
        self, source: dict[str, Any], rules: 'Rules', target: Any = None
    ) -> Any:
        """
        Applies all rules in a Rules object to transform the source dictionary into the target.

        Args:
            source (dict[str, Any]): The source data.
            rules (Rules): The set of rules to apply.
            target (Optional[Any], optional): The initial target data structure. Defaults to None.

        Returns:
            Any: The transformed target data structure.
        """
        if target is None:
            target = {} if isinstance(source, dict) else []
        for rule_name, rule in rules.rules.items():
            self.transform_dict(rule, source, target, self.mapping_dict)
        return target

    def transform(
        self,
        source_data: dict[str, Any],
        mapping_name: str | None = None,
        target_data: Any = None,
    ) -> Any:
        """
        Transforms the source data into the target data based on the specified mapping.

        Args:
            source_data (dict[str, Any]): The source JSON data.
            mapping_name (str): The name of the mapping to use. Default is None.
            target_data (Optional[Any], optional): The initial target data structure. Defaults to None.

        Raises:
            ValueError: If the specified mapping name does not exist.

        Returns:
            Any: The transformed target data structure.
        """
        if not mapping_name:
            mapping_name = list(source_data.keys())[0]
        try:
            if mapping_name not in self.mapping_dict:
                raise ValueError(
                    f"Mapping name '{mapping_name}' not found in the transformation dictionary"
                )

            mapping = self.mapping_dict[mapping_name]

            if target_data is None and any(
                rule.target.startswith('[') for rule in mapping.rules.values()
            ):
                target_data = []

            return self.dict_to_dict(source_data, mapping, target_data)
        except Exception as e:
            raise e
