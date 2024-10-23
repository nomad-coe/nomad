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

import {DateEditQuantity, DateTimeEditQuantity, TimeEditQuantity} from './DateTimeEditQuantity'
import { StringEditQuantity, URLEditQuantity } from './StringEditQuantity'
import {NumberEditQuantity} from './NumberEditQuantity'
import {EnumEditQuantity} from './EnumEditQuantity'
import {AutocompleteEditQuantity} from './AutocompleteEditQuantity'
import {BoolEditQuantity} from './BoolEditQuantity'
import {ActionEditQuantity} from './ActionEditQuantity'
import FileEditQuantity from './FileEditQuantity'
import RichTextEditQuantity from './RichTextEditQuantity'
import ReferenceEditQuantity from './ReferenceEditQuantity'
import AuthorEditQuantity from './AuthorEditQuantity'
import { RadioEnumEditQuantity } from './RadioEnumEditQuantity'
import QueryEditQuantity from "./QueryEditQuantity"
import {SliderEditQuantity} from "./SliderEditQuantity"

export const editQuantityComponents = {
  NumberEditQuantity: NumberEditQuantity,
  SliderEditQuantity: SliderEditQuantity,
  StringEditQuantity: StringEditQuantity,
  URLEditQuantity: URLEditQuantity,
  EnumEditQuantity: EnumEditQuantity,
  SelectEnumEditQuantity: EnumEditQuantity,
  RadioEnumEditQuantity: RadioEnumEditQuantity,
  AutocompleteEditQuantity: AutocompleteEditQuantity,
  BoolEditQuantity: BoolEditQuantity,
  ActionEditQuantity: ActionEditQuantity,
  FileEditQuantity: FileEditQuantity,
  DateTimeEditQuantity: DateTimeEditQuantity,
  DateEditQuantity: DateEditQuantity,
  TimeEditQuantity: TimeEditQuantity,
  RichTextEditQuantity: RichTextEditQuantity,
  ReferenceEditQuantity: ReferenceEditQuantity,
  AuthorEditQuantity: AuthorEditQuantity,
  QueryEditQuantity: QueryEditQuantity
}
