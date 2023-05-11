/* eslint-disable quotes */
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

import React, { memo, useEffect, useState } from 'react'
import PropTypes from 'prop-types'
import { PropertyCard } from './PropertyCard'
import { Box, Collapse, List, ListItem, ListItemText, Tooltip } from '@material-ui/core'
import { useApi } from '../../api'
import { useErrors } from '../../errors'
import Quantity from '../../Quantity'
import { Datatable, DatatableTable } from '../../datatable/Datatable'
import EntryDetails, { EntryRowActions } from '../EntryDetails'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ExpandLessIcon from '@material-ui/icons/ExpandLess'
import { Alert } from '@material-ui/lab'

const columns = [
  {
    key: 'entry_name',
    label: 'Name',
    align: 'left',
    sortable: false,
    render: entry => <Quantity quantity={'entry_name'} noWrap noLabel placeholder="unnamed" data={entry}/>
  },
  {key: 'entry_type', align: 'left', sortable: false, label: 'Type'},
  {
    key: 'mainfile',
    align: 'left',
    render: entry => <Quantity quantity={'mainfile'} noLabel noWrap withClipboard data={entry}/>
  },
  {
    key: 'upload_id',
    render: upload => <Quantity quantity={'upload_id'} noLabel noWrap withClipboard data={upload}/>
  }
]

const defaultColumns = [
  'entry_name',
  'entry_type',
  'mainfile'
]

const baseQuery = {
  'required': {
    'metadata': '*'
  },
  'owner': 'visible',
  'pagination': {
    'order_by': 'upload_create_time', 'order': 'desc', 'page_size': 20
  }
}

export const ReferenceTable = memo(({refPool, defaultColumns, noAction}) => {
  return <Datatable
    columns={columns} shownColumns={defaultColumns}
    data={refPool.map(entry => ({...entry.entry_metadata, ...entry}))}
  >
    <DatatableTable actions={!noAction && EntryRowActions} details={EntryDetails}>
    </DatatableTable>
  </Datatable>
})
ReferenceTable.propTypes = {
  refPool: PropTypes.arrayOf(
    PropTypes.shape({
      entry_metadata: PropTypes.object
    })
  ),
  defaultColumns: PropTypes.arrayOf(PropTypes.string),
  noAction: PropTypes.bool
}

const ReferenceUsingCard = memo(({index}) => {
  const [referencingPool, setReferencingPool] = useState(null)
  const [referencedPool, setReferencedPool] = useState(null)
  const [openReferencing, setOpenReferencing] = useState(true)
  const [openReferenced, setOpenReferenced] = useState(false)

  const {api} = useApi()
  const {raiseError} = useErrors()

  // the referencing list from the archive itself
  useEffect(() => {
    if (!index?.entry_references) {
      setReferencingPool(null)
      return
    }

    const targetEntryIds = new Set()
    for (const entry of index.entry_references) targetEntryIds.add(entry.target_entry_id)

    const referencing_query = {
      ...baseQuery,
      query: {
        'entry_id:any': [...targetEntryIds]
      }
    }
    api.query('entries', referencing_query, {noLoading: true})
      .then((data) => {
        setReferencingPool(data.data)
      })
      .catch(raiseError)
  }, [api, index, raiseError])

  const referencingList = referencingPool?.length > 0
    ? <ReferenceTable refPool={referencingPool} defaultColumns={defaultColumns}/>
    : <Alert severity="info">Not referencing other entries.</Alert>

  // load the referenced list of references from search
  useEffect(() => {
    const referenced_query = {
      ...baseQuery,
      ...{query: {'entry_references.target_entry_id': index.entry_id}}
    }
    api.query('entries', referenced_query, {noLoading: true}).then((data) => {
      setReferencedPool(data.data)
    }).catch(raiseError)
  }, [api, index, raiseError])

  // the referenced by list via search
  const referencedList = referencedPool?.length > 0
    ? <ReferenceTable refPool={referencedPool} defaultColumns={defaultColumns}/>
    : <Alert severity="info">Not referenced by other entries.</Alert>

  return (<PropertyCard title="Entry References">
    <List dense={true}>
      <ListItem button onClick={() => setOpenReferencing(!openReferencing)}>
        <Tooltip title="Click to expand/collapse">
          <ListItemText primaryTypographyProps={{variant: 'subtitle1'}}
                        primary="Referencing the following entries"/>
        </Tooltip>
        {openReferencing ? <ExpandLessIcon fontSize="small"/> : <ExpandMoreIcon fontSize="small"/>}
      </ListItem>
      <Collapse in={openReferencing} timeout="auto">
        <Box marginBottom={3}>
          {referencingList}
        </Box>
      </Collapse>
      <ListItem button onClick={() => setOpenReferenced(!openReferenced)}>
        <Tooltip
          title={(referencedPool?.length >= 20 ? 'Only showing the latest 20 entries. ' : '') + 'Click to expand/collapse'}>
          <ListItemText primaryTypographyProps={{variant: 'subtitle1'}}
                        primary="Referenced by the following entries"/>
        </Tooltip>
        {openReferenced ? <ExpandLessIcon fontSize="small"/> : <ExpandMoreIcon fontSize="small"/>}
      </ListItem>
      <Collapse in={openReferenced} timeout="auto">{referencedList}</Collapse>
    </List>
  </PropertyCard>)
})
ReferenceUsingCard.propTypes = {
  index: PropTypes.shape({
    entry_id: PropTypes.string.isRequired,
    entry_references: PropTypes.arrayOf(PropTypes.shape({
      target_entry_id: PropTypes.string.isRequired,
      target_path: PropTypes.string.isRequired
    }))
  }).isRequired
}

export default ReferenceUsingCard
