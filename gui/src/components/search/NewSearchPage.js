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
import { domainData } from '../domainData'

const help = `
This page allows you to **search** in NOMAD's data. The upper part of this page
gives you various options to enter and configure your search. The lower part
shows all data that fulfills your search criteria.

NOMAD's *domain-aware* search allows you to screen data by filtering based on
desired properties. This is different from basic *text-search* that traditional
search engines offer.

The search bar allows you to specify various quantity values that you want to
see in your results. This includes *authors*, *comments*, *atom labels*, *code name*,
*system type*, *crystal system*, *basis set types*, and *XC functionals*.

Alternatively, you can click *elements* or *metadata* to get a visual representation of
NOMAD's data as a periodic table or metadata charts. You can click the various
visualization elements to filter for respective quantities.

The visual representations show metrics for all data that fit your criteria.
You can display *entries* (i.e. code runs), *unique entries*, and *datasets*.
Other more specific metrics might be available.

Some quantities have no autocompletion for their values. You can still search for them,
if you know exactly what you are looking for. To search for a particular entry by its id
for example, type \`calc_id=<the_id>\` and press entry (or select the respective item from the menu).
The usable *hidden* quantities are: ${Object.keys(domainData.dft.additionalSearchKeys).map(key => `\`${key}\``).join(', ')}.

The results tabs gives you a quick overview of all entries and datasets that fit your search.
You can click entries to see more details, download data, see the archive, etc. The *entries*
tab displays individual entries (i.e. code runs), the *grouped entries* tab will group
entries with similar metadata (it will group entries for the same material from the
  same user). The *dataset* tab, shows entry curated by user created datasets. You can
  click on datasets for a search page that will only display entries from the respective
  dataset.

The table columns can be configured. The *entries* tab also supports sorting. Selected
entries (or all entries) can be downloaded. The download will contain all user provided
raw calculation input and output files.

You can click entries to see more details about them. The details button will navigate
you to an entry's page. This entry page will show more metadata, raw files, the
entry's archive, and processing logs.
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
