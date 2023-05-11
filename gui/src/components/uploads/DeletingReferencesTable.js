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
import React, {useMemo} from 'react'
import PropTypes from 'prop-types'
import {DialogContentText} from '@material-ui/core'
import {Datatable} from '../datatable/Datatable'
import Quantity from '../Quantity'

const columns = [
  {key: 'mainfile', align: 'left', label: 'Mainfile'},
  {key: 'upload_id', align: 'left', label: 'Upload Id', render: data => <Quantity quantity={'upload_id'} noLabel noWrap withClipboard data={data}/>}
]

export const DeletingReferencesTable = React.memo(({entryReferences, brokenEntries, message}) => {
  const normalizedEntryReferences = useMemo(() => {
    return entryReferences.map(reference => ({...reference, mainfile: reference.target_mainfile, upload_id: reference.target_upload_id}))
  }, [entryReferences])

  return (
    brokenEntries.length > 0 ? <React.Fragment>
    <DialogContentText>
      <span style={{color: 'red'}}>{message || 'Warning: There are some entries which are referenced somewhere else.'}</span>
    </DialogContentText>
    <DialogContentText>
      References Being Deleted
      <Datatable columns={columns} data={normalizedEntryReferences}/>
    </DialogContentText>
    <DialogContentText>
      Broken Entries
      <Datatable columns={columns} data={brokenEntries}/>
    </DialogContentText>
  </React.Fragment> : null)
})
DeletingReferencesTable.propTypes = {
  entryReferences: PropTypes.arrayOf(PropTypes.object).isRequired,
  brokenEntries: PropTypes.arrayOf(PropTypes.object).isRequired,
  message: PropTypes.string
}

export default DeletingReferencesTable
