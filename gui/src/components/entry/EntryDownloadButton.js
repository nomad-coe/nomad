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
import React, { useState } from 'react'
import PropTypes from 'prop-types'
import { useErrors } from '../errors'
import { apiBase } from '../../config'
import { Tooltip, IconButton, Menu, MenuItem } from '@material-ui/core'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import { useApi } from '../api'
import { toAPIFilter } from '../search/SearchContext'

const EntryDownloadButton = React.memo(function EntryDownloadButton(props) {
  const {tooltip, disabled, buttonProps, dark, query} = props
  const {api, user} = useApi()
  const {raiseError} = useErrors()

  const [preparingDownload, setPreparingDownload] = useState(false)
  const [anchorEl, setAnchorEl] = useState(null)

  const download = (choice) => {
    let queryStringData = toAPIFilter(query)
    const owner = query.visibility || 'visible'
    const openDownload = (token) => {
      const url = `${apiBase}/v1/entries/${choice}/download?owner=${owner}&signature_token=${token}&json_query=${JSON.stringify(queryStringData)}`
      window.location.assign(url)
    }

    if (user) {
      setPreparingDownload(true)
      api.get('/auth/signature_token')
        .then(response => {
          const token = response.signature_token
          openDownload(token)
        })
        .catch(raiseError)
        .finally(setPreparingDownload(false))
    } else {
      openDownload()
    }
  }

  const handleClick = event => {
    event.stopPropagation()
    setAnchorEl(event.currentTarget)
  }

  const handleSelect = (choice) => {
    setAnchorEl(null)
    download(choice)
  }

  const handleClose = () => {
    setAnchorEl(null)
  }

  return <React.Fragment>
    <IconButton
      {...buttonProps}
      disabled={disabled || preparingDownload}
      onClick={handleClick}
      style={dark ? {color: 'white'} : null}
    >
      <Tooltip title={tooltip || 'Download'}>
        <DownloadIcon />
      </Tooltip>
    </IconButton>
    <Menu
      anchorEl={anchorEl}
      open={Boolean(anchorEl)}
      onClose={handleClose}
    >
      <MenuItem onClick={() => handleSelect('raw')}>Raw uploaded files</MenuItem>
      <MenuItem onClick={() => handleSelect('archive')}>Processed data</MenuItem>
    </Menu>
  </React.Fragment>
})
EntryDownloadButton.propTypes = {
  /**
   * The query that defines what to download.
   */
  query: PropTypes.object.isRequired,
  /**
   * A tooltip for the button
   */
  tooltip: PropTypes.string,
  /**
   * Whether the button is disabled
   */
  disabled: PropTypes.bool,
  /**
   * Properties forwarded to the button.
   */
  buttonProps: PropTypes.object,
  dark: PropTypes.bool
}

export default EntryDownloadButton
