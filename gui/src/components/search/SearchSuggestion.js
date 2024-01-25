/*
 * Copyright The NOMAD Authors.
 *
 * This file is part of NOMAD. See https://nomad-lab.eu for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * Represents a suggestable option in the search bar. In TypeScript this would
 * be an interface.
 */
export class SearchSuggestion {
  input
  type
  history

  constructor(param) {
    this.input = param.input
    this.type = param.type
    this.history = param.history
    this.key = `${this.input}-${this.type}`
  }
}

/**
 * Enum for types of suggestions. In TypeScript this would be an Enum.
 */
export const SuggestionType = {
  Freetext: 'free_text',
  Equality: 'equality',
  RangeHalfBounded: 'range_half_bounded',
  RangeBounded: 'range_bounded',
  Existence: 'existence',
  Name: 'name'
}
