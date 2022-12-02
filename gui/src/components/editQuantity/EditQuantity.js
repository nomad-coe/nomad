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

import {DateTimeEditQuantity} from '../editQuantity/DateTimeEditQuantity'
import { StringEditQuantity, URLEditQuantity } from '../editQuantity/StringEditQuantity'
import {NumberEditQuantity} from '../editQuantity/NumberEditQuantity'
import {EnumEditQuantity} from '../editQuantity/EnumEditQuantity'
import {AutocompleteEditQuantity} from '../editQuantity/AutocompleteEditQuantity'
import {BoolEditQuantity} from '../editQuantity/BoolEditQuantity'
import FileEditQuantity from '../editQuantity/FileEditQuantity'
import RichTextEditQuantity from '../editQuantity/RichTextEditQuantity'
import ReferenceEditQuantity from '../editQuantity/ReferenceEditQuantity'
import AuthorEditQuantity from '../editQuantity/AuthorEditQuantity'
import { RadioEnumEditQuantity } from '../editQuantity/RadioEnumEditQuantity'

export const editQuantityComponents = {
  NumberEditQuantity: NumberEditQuantity,
  StringEditQuantity: StringEditQuantity,
  URLEditQuantity: URLEditQuantity,
  EnumEditQuantity: EnumEditQuantity,
  SelectEnumEditQuantity: EnumEditQuantity,
  RadioEnumEditQuantity: RadioEnumEditQuantity,
  AutocompleteEditQuantity: AutocompleteEditQuantity,
  BoolEditQuantity: BoolEditQuantity,
  FileEditQuantity: FileEditQuantity,
  DateTimeEditQuantity: DateTimeEditQuantity,
  RichTextEditQuantity: RichTextEditQuantity,
  ReferenceEditQuantity: ReferenceEditQuantity,
  AuthorEditQuantity: AuthorEditQuantity
}
