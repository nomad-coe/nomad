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
import FileSaver from 'file-saver'
import { useErrors } from './errors'
import { apiBase } from '../config'
import { Tooltip, IconButton, Menu, MenuItem } from '@material-ui/core'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import { useApi } from './api'

const DownloadButton = React.memo(function DownloadButton(props) {
  const {tooltip, disabled, buttonProps, dark, query} = props
  const {api, user} = useApi()
  const {raiseError} = useErrors()

  const [preparingDownload, setPreparingDownload] = useState(false)
  const [anchorEl, setAnchorEl] = useState(null)

  const handleClick = event => {
    event.stopPropagation()
    setAnchorEl(event.currentTarget)
  }

  const handleSelect = (choice) => {
    setAnchorEl(null)

    const urlSearchParams = new URLSearchParams()
    Object.keys(query).forEach(key => {
      const value = query[key]
      if (Array.isArray(value)) {
        value.forEach(item => urlSearchParams.append(key, item))
      } else {
        urlSearchParams.append(key, value)
      }
    })

    const openDownload = () => {
      const url = `${apiBase}/${choice}/${choice === 'archive' ? 'download' : 'query'}?${urlSearchParams.toString()}`
      FileSaver.saveAs(url, `nomad-${choice}-download.zip`)
    }

    if (user) {
      setPreparingDownload(true)
      api.get('/auth/signature_token')
        .then(response => {
          urlSearchParams.append('signature_token', response.signature_token)
          openDownload()
        })
        .catch(raiseError)
        .finally(setPreparingDownload(false))
    } else {
      openDownload()
    }
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
      <MenuItem onClick={() => handleSelect('archive')}>NOMAD Archive files</MenuItem>
    </Menu>
  </React.Fragment>
})
DownloadButton.propTypes = {
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

export default DownloadButton
