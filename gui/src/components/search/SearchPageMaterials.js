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
import React from 'react'
import Search from './Search'
import { SearchContext } from './SearchContext'

export const help = `
This page allows you to **search materials** within NOMAD. NOMAD can
automatically detect the material from individual entries and can then group the
data by using these detected materials. This allows you to search individual
materials which have properties that are aggregated from several entries.

The search page consists of three main elements: the filter panel, the search
bar, and the result list.

The filter panel on the left allows you to graphically explore and enter
different search filters. It also gives a visual indication of the currently
active search filters for each category. This is a good place to start exploring
the available search filters and their meaning.

The search bar allows you to specify filters by typing them in and pressing
enter. You can also start by simply typing keywords of interest, which will
toggle a list of suggestions. For numerical data you can also use range queries,
e.g. \`0.0 < band_gap <= 0.1\`.

The units used in the filter panel and in the queries can be changed
using the **units** button on the top right corner. When using the search bar,
you can also specify a unit by typing the unit abbreviations, e.g. \`band_gap >=
0.1 Ha\`.

Notice that by default the properties that you search can be combined from
several different entries. If instead you wish to search for a material with an
individual entry fullfilling your search criteria, uncheck the **combine results
from several entries**-checkbox.

The result list on the right is automatically updated according to the filters
you have specified. You can scroll through the available items and load more
results as you go. Here you can also change the sorting of the results, modify
the displayed columns and access individual materials. The ellipsis button shown
for each material will navigate you into the material overview page within the
NOMAD Encyclopedia. This page will show a more detailed overview for that
specific material.
`

const SearchPageMaterials = React.memo(() => {
  return <SearchContext resource="materials">
    <Search/>
  </SearchContext>
})

export default SearchPageMaterials
