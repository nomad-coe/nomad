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
import React, { useContext, useState, useEffect, useMemo } from 'react'
import PropTypes from 'prop-types'
import { cloneDeep } from 'lodash'
import { Typography, makeStyles } from '@material-ui/core'
import { ui } from '../../config'
import { errorContext } from '../errors'
import { useApi } from '../api'
import SearchPage from '../search/SearchPage'
import { SearchContext } from '../search/SearchContext'
import { DOI } from './DOI'

export const help = `
This page allows you to **inspect** and **download** NOMAD datasets. It also allows you
to explore a dataset with similar controls that the search page offers.
`

// Use the same context as in the global entries search, but with free-text
// query enabled
const context = cloneDeep(ui?.apps?.options?.entries)
context.search_syntaxes.exclude = undefined

const useStyles = makeStyles(theme => ({
  header: {
    display: 'flex',
    flexDirection: 'column'
  }
}))
const DatasetPage = React.memo(({match}) => {
  const styles = useStyles()
  const [dataset, setDataset] = useState()
  const {raiseError} = useContext(errorContext)
  const {api} = useApi()

  // Router provides the URL parameters via props, here we read the dataset ID.
  const datasetId = match?.params?.datasetId
  const datasetFilter = useMemo(() => ({'datasets.dataset_id': datasetId}), [datasetId])

  // Fetch the dataset information from API.
  useEffect(() => {
    api.datasets(datasetId)
      .then(setDataset)
      .catch(error => {
        setDataset(undefined)
        raiseError(error)
      })
  }, [datasetId, api, raiseError])

  // Show a customized search page for this dataset. Basic dataset information
  // shown in the header.
  return dataset
    ? <SearchContext
      resource={context?.resource}
      initialPagination={context?.pagination}
      initialColumns={context?.columns}
      initialRows={context?.rows}
      initialFilterMenus={context?.filter_menus}
      initialFilters={context?.filters}
      initialFiltersLocked={datasetFilter}
      initialDashboard={context?.dashboard}
      initialSearchSyntaxes={context?.search_syntaxes}
      id='dataset-search'
    >
      <SearchPage header={
        <div className={styles.header}>
          <Typography variant="h4">
            {dataset.dataset_name || (dataset.isEmpty && 'Empty or non existing dataset') || 'loading ...'}
          </Typography>
          <Typography>
            dataset{dataset.doi ? <span>, with DOI <DOI doi={dataset.doi} /></span> : ''}
          </Typography>
        </div>}
      />
    </SearchContext>
    : null
})
DatasetPage.propTypes = {
  match: PropTypes.object
}

export default DatasetPage
