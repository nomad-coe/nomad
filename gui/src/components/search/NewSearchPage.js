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
import React, { useContext } from 'react'
import { apiContext } from '../api'
import NewSearch from './NewSearch'

const help = `
This page allows you to **search** in NOMAD's data. NOMAD's *domain-aware*
search allows you to screen data by filtering based on desired properties. This
is different from basic *text-search* that traditional search engines offer.

The search page consists of three main elements: the filter panel on the left,
the search bar, and the result list on the right.

The filter panel on the left allows you to graphically explore and enter
different search filters. It also gives a visual indication of the currently
active search filters for each category. This is a good place to start exploring
the available search filters and their meaning.

The search bar allows you to specify filters by typing them in and pressing
enter. You can also start by simply typing keywords of interest, which will
toggle a list of possible suggestions. All of the filters that are available
through the left panel are also available in this search bar. For numerical data
you can also use range queries, e.g. \`0.0 < band_gap <= 0.1\`. The units used in the
queries can be changed in the settings.

The result list is automatically updated according to the filters you have
specified. You can browse through the results by simply scrolling through the
available items. Here you can also change the sorting of the results, modify the
displayed columns, access individual entries or even download selections of the
data. The results tabs gives you a quick overview of all entries and datasets
that fit your search and it is automatically updated based on your filters. You
can browse through all of the results by scrolling down the list. Here you can
also change the sorting of the results, modify the displayed columns, access
individual entries or even download selections of the data. The arrow button
shown for each entry will navigate you to that entry's page.  This entry page
will show more metadata, raw files, the entry's archive, and processing logs.
`
export {help}

export default function NewSearchPage() {
  const {user} = useContext(apiContext)
  const withoutLogin = ['all', 'public']

  return <NewSearch
    initialOwner="public"
    ownerTypes={['public', 'visible'].filter(key => user || withoutLogin.indexOf(key) !== -1)}
    showDisclaimer
  />
}
