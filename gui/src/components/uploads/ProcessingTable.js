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

import React, { useMemo, useEffect } from 'react'
import PropTypes from 'prop-types'

import EntryList, { EntryListUnstyled } from '../search/EntryList'
import { Paper } from '@material-ui/core'
import { authorList, nameList } from '../../utils'
import WithButton from '../utils/WithButton'

const defaultSelectedColumns = ['mainfile', 'parser', 'process_status', 'complete_time']

export default function ProcessingTable({data, onPaginationChange}) {
  const columns = useMemo(() => {
    const otherColumns = {...EntryListUnstyled.defaultColumns}
    Object.keys(otherColumns).forEach(key => {
      otherColumns[key] = {
        ...otherColumns[key],
        supportsSort: false
      }
    })

    return {
      entry_id: {
        label: 'Id',
        description: 'The unique NOMAD id for this entry.',
        render: entry => <WithButton clipboard={entry.entry_id}>
          {entry.entry_id}
        </WithButton>,
        tableCellStyle: {maxWidth: 300}
      },
      mainfile: {
        ...EntryListUnstyled.defaultColumns.mainfile,
        supportsSort: true
      },
      parser: {
        label: 'Parser',
        supportsSort: true,
        description: 'The parser that was used to process this entry.',
        render: entry => entry.parser.replace('parsers/', '')
      },
      process_status: {
        label: 'Processing',
        supportsSort: true,
        description: 'Details on the processing of this entry.',
        render: entry => `${entry.process_status}`
      },
      complete_time: {
        label: 'Last processing',
        description: 'The last time this entry has completed its processing.',
        render: entry => new Date(entry.complete_time).toLocaleString()
      },
      comment: EntryListUnstyled.defaultColumns.comment,
      references: {
        ...EntryListUnstyled.defaultColumns.references,
        supportsSort: false
      },
      authors: {
        label: 'Authors',
        render: entry => authorList(entry),
        description: 'The uploader and co-authors of this entry.'
      },
      owners: {
        label: 'Owner',
        render: entry => nameList(entry.owners || []),
        description: 'The uploader and shared with users.'
      },
      datasets: EntryListUnstyled.defaultColumns.datasets
    }
  }, [])

  const handleChange = (changes) => {
    changes.page_size = changes.page_size || changes.per_page || data.pagination.page_size
    if (changes.order !== undefined) {
      changes.order = changes.order === 1 ? 'asc' : 'desc'
    }
    onPaginationChange({
      page_size: data.pagination.page_size,
      page: data.pagination.page,
      order_by: data.pagination.order_by,
      order: data.pagination.order,
      ...changes
    })
  }

  useEffect(() => {
    data.data.forEach(entry => {
      Object.assign(entry, entry.entry_metadata || {
        references: [], authors: [], owners: [], datasets: []
      })
    })
  }, [data.data])

  const editable = !data.upload.process_running

  return <Paper>
    <EntryList
      title={`${data.pagination.total} entries processed`}
      query={{upload_id: [data.upload.upload_id]}}
      columns={columns}
      selectedColumns={defaultSelectedColumns}
      selectedColumnsKey={null}
      editable={editable}
      onEdit={() => handleChange({})}
      editUserMetadataDialogProps={{withoutLiftEmbargo: !data.upload.published}}
      data={data}
      onChange={handleChange}
      // actions={actions}
      // showEntryActions={entry => entry.processed || !running}
      showEntryActions={entry => !entry.process_running}
      entryPagePathPrefix={`/uploads/${data.upload.upload_id}`}
      per_page={data.pagination.page_size}
      page={data.pagination.page}
      order={data.pagination.order === 'asc' ? 1 : 0}
      order_by={data.pagination.order_by}
      // TODO this is a hack, gave up trying to understand CSS width in table cells
      entryDetailsWidth={903}
    />
  </Paper>
}
ProcessingTable.propTypes = {
  data: PropTypes.object,
  onPaginationChange: PropTypes.func
}
