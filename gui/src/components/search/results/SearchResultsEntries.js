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
import React, { useState, useMemo } from 'react'
import PropTypes from 'prop-types'
import { Link } from '@material-ui/core'
import EntryDownloadButton from '../../entry/EntryDownloadButton'
import EntryDetails, { EntryRowActions, Published } from '../../entry/EntryDetails'
import { authorList, pluralize, formatInteger, formatTimestamp } from '../../../utils'
import {
  addColumnDefaults,
  Datatable, DatatableLoadMorePagination, DatatableTable,
  DatatableToolbar, DatatableToolbarActions } from '../../datatable/Datatable'
import { useSearchContext } from '../SearchContext'

const columns = [
  {key: 'results.material.chemical_formula_hill', label: 'Formula', align: 'left'},
  {key: 'results.method.method_name'},
  {key: 'results.method.simulation.program_name'},
  {key: 'results.method.simulation.dft.basis_set_name'},
  {key: 'results.method.simulation.dft.xc_functional_type', label: 'XC functional type'},
  {key: 'results.material.structural_type'},
  {key: 'results.material.symmetry.crystal_system'},
  {key: 'results.material.symmetry.space_group_symbol'},
  {key: 'results.material.symmetry.space_group_number'},
  {key: 'mainfile', align: 'left'},
  {
    key: 'upload_create_time',
    label: 'Upload time',
    align: 'left',
    render: row => row?.upload_create_time
      ? formatTimestamp(row.upload_create_time)
      : <i>no upload time</i>
  },
  {key: 'authors', render: row => authorList(row), align: 'left', sortable: false},
  {key: 'comment', sortable: false, align: 'left'},
  {
    key: 'references',
    sortable: false,
    align: 'left',
    render: row => {
      const refs = row.references || []
      if (refs.length > 0) {
        return (
          <div style={{display: 'inline'}}>
            {refs.map((ref, i) => <span key={ref}>
              <Link href={ref}>{ref}</Link>{(i + 1) < refs.length ? ', ' : <React.Fragment/>}
            </span>)}
          </div>
        )
      } else {
        return <i>no references</i>
      }
    }
  },
  {
    key: 'datasets',
    align: 'left',
    sortable: false,
    render: entry => {
      const datasets = entry.datasets || []
      if (datasets.length > 0) {
        return datasets.map(dataset => dataset.dataset_name).join(', ')
      } else {
        return <i>no datasets</i>
      }
    }
  },
  {
    key: 'published',
    label: 'Access',
    render: (entry) => <Published entry={entry} />
  }
]

addColumnDefaults(columns)

const defaultSelectedColumns = [
  'results.material.chemical_formula_hill',
  'results.method.method_name',
  'upload_create_time',
  'authors'
]

const SearchResultsEntries = React.memo((props) => {
  const { useQuery } = useSearchContext()
  const [selected, setSelected] = useState([])
  const searchQuery = useQuery()
  const {pagination, data} = props

  const query = useMemo(() => {
    if (selected === 'all') {
      return searchQuery
    }

    return {entry_id: selected.map(data => data.entry_id)}
  }, [selected, searchQuery])

  return <Datatable
    columns={columns} shownColumns={defaultSelectedColumns} {...props}
    selected={selected} onSelectedChanged={setSelected}
  >
    <DatatableToolbar title={`${formatInteger(data.length)}/${pluralize('result', pagination.total, true, true, 'search')}`}>
      <DatatableToolbarActions selection>
        <EntryDownloadButton tooltip="Download files" query={query} />
      </DatatableToolbarActions>
    </DatatableToolbar>
    <DatatableTable actions={EntryRowActions} details={EntryDetails}>
      <DatatableLoadMorePagination color="primary">load more</DatatableLoadMorePagination>
    </DatatableTable>
  </Datatable>
})
SearchResultsEntries.propTypes = {
  pagination: PropTypes.object,
  data: PropTypes.arrayOf(PropTypes.object)
}

export default SearchResultsEntries
