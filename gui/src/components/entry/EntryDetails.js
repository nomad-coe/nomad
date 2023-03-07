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
import PropTypes from 'prop-types'
import Quantity from '../Quantity'
import { Formula } from './properties/MaterialCard'
import { Tooltip, IconButton } from '@material-ui/core'
import { useApi } from '../api'
import { EntryButton } from '../nav/Routes'
import DetailsIcon from '@material-ui/icons/MoreHoriz'
import PublicIcon from '@material-ui/icons/Public'
import UploaderIcon from '@material-ui/icons/AccountCircle'
import SharedIcon from '@material-ui/icons/SupervisedUserCircle'
import PrivateIcon from '@material-ui/icons/VisibilityOff'
import { makeStyles } from '@material-ui/core/styles'

export const MethodMetadata = React.memo(({data}) => {
  const methodQuantities = []
  const addMethodQuantities = (obj, parentKey, hidden = []) => {
    const children = {}
    Object.keys(obj).forEach(key => {
      if (hidden.includes(key)) {
        return
      }
      const value = obj[key]
      if (Array.isArray(value) || typeof value === 'string') {
        if (value.length > 0) {
          methodQuantities.push({quantity: `${parentKey}.${key}`})
        }
      } else if (value instanceof Object) {
        children[key] = value
      }
    })
    Object.keys(children).forEach(key => addMethodQuantities(children[key], `${parentKey}.${key}`))
  }
  if (data?.results?.method) {
    addMethodQuantities(data.results.method, 'results.method')
  }
  if (data?.results?.eln) {
    addMethodQuantities(data.results.eln, 'results.eln', ['names', 'sections', 'descriptions'])
  }

  methodQuantities.push({quantity: 'entry_type', label: 'type'})
  methodQuantities.push({quantity: 'entry_name', label: 'name'})

  return <Quantity flex>
    {methodQuantities.map(({...quantityProps}) => (
      <Quantity
        key={quantityProps.quantity}
        {...quantityProps}
        noWrap
        data={data}
        hideIfUnavailable
      />
    ))}
  </Quantity>
})
MethodMetadata.propTypes = {
  data: PropTypes.object
}

export const DomainMetadata = React.memo(({data}) => {
  return <>
    <Quantity flex>
      <Formula data={data} />
    </Quantity>
    <MethodMetadata data={data} />
  </>
})
DomainMetadata.propTypes = {
  data: PropTypes.object
}

export const UserMetadata = React.memo(({data}) => {
  return (
    <div>
      <Quantity quantity='comment' data={data} />
      <Quantity quantity='references' data={data}/>
      <Quantity quantity='authors' data={data}/>
      <Quantity quantity='datasets' data={data}/>
    </div>
  )
})
UserMetadata.propTypes = {
  data: PropTypes.object.isRequired
}

export const EntryIds = React.memo(({data}) => {
  return (
    <div>
      <Quantity column >
        {/* <Quantity quantity="pid" label='PID' placeholder="not yet assigned" noWrap data={data} withClipboard /> */}
        <Quantity quantity="entry_id" noWrap withClipboard data={data} />
        <Quantity quantity="raw_id" noWrap withClipboard data={data} />
        <Quantity quantity="external_id" noWrap withClipboard data={data} />
        <Quantity quantity="mainfile" noWrap ellipsisFront data={data} withClipboard />
        <Quantity quantity="upload_id" data={data}/>
      </Quantity>
    </div>
  )
})
EntryIds.propTypes = {
  data: PropTypes.object.isRequired
}

export function Published(props) {
  const {user} = useApi()
  const {entry} = props
  if (entry.published) {
    if (entry.with_embargo) {
      if (user && entry.main_author.user_id === user.sub) {
        if (entry.viewers.length === 1) {
          return <Tooltip title="published with embargo by you and only accessible by you">
            <UploaderIcon color="error" />
          </Tooltip>
        } else {
          return <Tooltip title="published with embargo by you and only accessible to you and the specified coauthors and reviewers">
            <SharedIcon color="error" />
          </Tooltip>
        }
      } else if (user && entry.coauthors.find(user => user.user_id === user.sub)) {
        return <Tooltip title="published with embargo and visible to you as a coauthor">
          <SharedIcon color="error" />
        </Tooltip>
      } else if (user && entry.reviewers.find(user => user.user_id === user.sub)) {
        return <Tooltip title="published with embargo and visible to you as a reviewer">
          <SharedIcon color="error" />
        </Tooltip>
      } else {
        if (user) {
          return <Tooltip title="published with embargo and not accessible by you">
            <PrivateIcon color="error" />
          </Tooltip>
        } else {
          return <Tooltip title="published with embargo and might become accessible after login">
            <PrivateIcon color="error" />
          </Tooltip>
        }
      }
    } else {
      return <Tooltip title="published and accessible by everyone">
        <PublicIcon color="primary" />
      </Tooltip>
    }
  } else {
    return <Tooltip title="you have not published this entry yet">
      <UploaderIcon color="error"/>
    </Tooltip>
  }
}
Published.propTypes = {
  entry: PropTypes.object.isRequired
}

export const VisitEntryAction = React.memo(function VisitEntryAction({data, ...props}) {
  const {user} = useApi()
  const hide = (data.with_embargo && !user && !data.viewers.find(viewer => viewer.user_id === user.sub)) || data.process_running
  if (hide) {
    return null
  }

  // The portal is disabled for this tooltip because this button causes a
  // navigation that otherwise leaves the popup opened (the Tooltip state does
  // not get updated since the page is cached and a new page is shown
  // immediately).
  return <Tooltip PopperProps={{disablePortal: true}} title="Go to the entry page">
    <EntryButton
      {...props}
      entryId={data.entry_id}
    />
  </Tooltip>
})
VisitEntryAction.propTypes = {
  data: PropTypes.object.isRequired
}

export const EntryRowActions = React.memo((props) => {
  return <VisitEntryAction {...props} component={IconButton}><DetailsIcon/></VisitEntryAction>
})

const useEntryDetailsStyles = makeStyles(theme => ({
  entryDetails: {
    paddingTop: theme.spacing(2),
    paddingLeft: theme.spacing(4),
    paddingRight: theme.spacing(4)
  },
  entryDetailsContents: {
    display: 'flex',
    width: '100%',
    margin: '0'
  },
  entryDetailsRow: {
    paddingRight: theme.spacing(3)
  },
  entryDetailsActions: {
    display: 'flex',
    flexBasis: 'auto',
    flexGrow: 0,
    flexShrink: 0,
    justifyContent: 'flex-end',
    marginBottom: theme.spacing(1),
    marginTop: theme.spacing(2)
  }
}))

export const EntryDetails = React.memo(({data}) => {
  const classes = useEntryDetailsStyles()

  return (
    <div className={classes.entryDetails}>
      <div className={classes.entryDetailsContents}>
        <div className={classes.entryDetailsRow}>
          <DomainMetadata data={data} />
        </div>

        <div className={classes.entryDetailsRow} style={{flexGrow: 1, minWidth: 'fit-content'}}>
          <Quantity className={classes.entryDetailsRow} column>
            <Quantity quantity='comment' data={data} />
            <Quantity quantity='references' data={data} />
            <Quantity quantity='authors' data={data} />
            <Quantity quantity='datasets' data={data} />
          </Quantity>
        </div>

        <div className={classes.entryDetailsRow} style={{maxWidth: '33%', paddingRight: 0}}>
          <Quantity column >
            <Quantity quantity="entry_id" data={data} />
            <Quantity quantity="mainfile" data={data} />
            <Quantity quantity="upload_id" data={data} />
            <Quantity quantity="raw_id" data={data} />
            <Quantity quantity="external_id" data={data} />
          </Quantity>
        </div>
      </div>

      <div className={classes.entryDetailsActions}>
        <VisitEntryAction color="primary" data={data}>
          Go to the entry page
        </VisitEntryAction>
      </div>
    </div>
  )
})
EntryDetails.propTypes = {
  data: PropTypes.object.isRequired
}

export default EntryDetails
