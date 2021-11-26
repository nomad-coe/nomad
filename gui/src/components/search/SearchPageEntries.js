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
This page allows you to **search entries** within NOMAD. Entries represent
individual calculations or experiments that have bee uploaded into NOMAD.

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

Notice that the units used in the filter panel and in the queries can be changed
using the **units** button on the top right corner. When using the search bar,
you can also specify a unit by typing the unit abbreviations, e.g. \`band_gap >=
0.1 Ha\`

The result list on the right is automatically updated according to the filters
you have specified. You can browse through the results by scrolling through the
available items and loading more results as you go. Here you can also change the
sorting of the results, modify the displayed columns, access individual entries
or even download the raw data or the archive document by selecting individual
entries and pressing the download button that appears. The ellipsis button shown
for each entry will navigate you to that entry's page. This entry page will show
more metadata, raw files, the entry's archive, and processing logs.
`

const SearchPageEntries = React.memo(() => {
  return <SearchContext resource="entries">
    <Search/>
  </SearchContext>
})

export default SearchPageEntries
