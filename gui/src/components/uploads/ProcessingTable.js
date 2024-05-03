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

import React, {useMemo, useState} from 'react'
import PropTypes from 'prop-types'
import { Paper, Link } from '@material-ui/core'
import EntryDetails, { EntryRowActions } from '../entry/EntryDetails'
import {
  addColumnDefaults,
  Datatable, DatatablePagePagination, DatatableTable,
  DatatableToolbar, DatatableToolbarActions } from '../datatable/Datatable'
import EntryDownloadButton from '../entry/EntryDownloadButton'
import DeleteEntriesButton from './DeleteEntriesButton'
import Quantity from '../Quantity'
import EditMetaDataDialog from './EditMetaDataDialog'
import {authorList, formatTimestamp, pluralize} from '../../utils'
import { useUploadPageContext } from './UploadPageContext'
import { defaultFilterData } from '../search/FilterRegistry'

const columns = [
  {
    key: 'entry_name',
    label: 'Name',
    align: 'left',
    sortable: false,
    render: entry => <Quantity quantity={'entry_name'} noWrap noLabel placeholder="unnamed" data={entry} />
  },
  {key: 'entry_type', align: 'left', sortable: false, label: 'Type'},
  {
    key: 'mainfile',
    align: 'left',
    render: entry => <Quantity quantity={'mainfile'} noLabel noWrap withClipboard data={entry}/>
  },
  {
    key: 'entry_id',
    align: 'left',
    sortable: false,
    render: entry => <Quantity quantity={'entry_id'} noLabel noWrap withClipboard data={entry}/>
  },
  {key: 'parser_name', align: 'left'},
  {key: 'authors', align: 'left', render: row => authorList(row)},
  {key: 'process_status', align: 'left'},
  {
    label: 'Modified',
    key: 'complete_time',
    align: 'left',
    sortable: false,
    render: entry => formatTimestamp(entry.complete_time)
  },
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
  }
]

addColumnDefaults(columns, defaultFilterData)

const defaultSelectedColumns = [
  'entry_name',
  'entry_type',
  'mainfile',
  'process_status'
]

export default function ProcessingTable(props) {
  const [selected, setSelected] = useState(new Set())
  const {pagination, customTitle} = props
  const {upload, isEditable} = useUploadPageContext()

  const selectedQuery = useMemo(() => {
    if (selected === 'all') {
      return {'upload_id': upload.upload_id}
    }

    return {entry_id: [...selected]}
  }, [selected, upload])

  return <Paper data-testid={'processing-table'}>
    <Datatable
      columns={columns} shownColumns={defaultSelectedColumns} {...props}
      selected={selected} getId={option => option.entry_id} onSelectedChanged={setSelected}
    >
      <DatatableToolbar title={pluralize((customTitle || 'search result'), pagination.total, true)}>
        <DatatableToolbarActions selection>
          <EntryDownloadButton tooltip="Download files" query={selectedQuery} />
          {isEditable &&
            <EditMetaDataDialog
              isIcon
              selectedEntries={selectedQuery}
            />}
          {isEditable &&
            <DeleteEntriesButton
              isIcon
              selectedEntries={selectedQuery}
              selectedCount={selected === 'all' ? pagination.total : selected?.size}
              setSelected={setSelected}
            />}
        </DatatableToolbarActions>
      </DatatableToolbar>
      <DatatableTable actions={EntryRowActions} details={EntryDetails}>
        <DatatablePagePagination pageSizeValues={[5, 10, 50, 100]} />
      </DatatableTable>
    </Datatable>
  </Paper>
}
ProcessingTable.propTypes = {
  data: PropTypes.arrayOf(PropTypes.object).isRequired,
  pagination: PropTypes.object.isRequired,
  onPaginationChanged: PropTypes.func.isRequired,
  customTitle: PropTypes.string
}
