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
import { useErrors } from '../errors'
import { apiBase } from '../../config'
import { makeStyles, Tooltip } from '@material-ui/core'
import { useApi } from '../api'

const useStyles = makeStyles(theme => ({
  root: {}
}))

export const Download = React.memo(function Download(props) {
  const classes = useStyles(props)
  const {url, fileName, component, children, disabled, color, size, tooltip} = props
  const {api, user} = useApi()
  const {raiseError} = useErrors()

  const [preparingDownload, setPreparingDownload] = useState(false)

  const handleClick = () => {
    let fullUrl = `${apiBase}/v1/${url}`
    let downloadUrl = fullUrl
    if (user) {
      setPreparingDownload(true)
      api.get('/auth/signature_token')
        .then(response => {
          if (fullUrl.startsWith('/')) {
            fullUrl = `${window.location.origin}${fullUrl}`
          }
          const downloadUrl = new URL(fullUrl)
          downloadUrl.searchParams.append('signature_token', response.signature_token)
          FileSaver.saveAs(downloadUrl.href, fileName)
        })
        .catch(raiseError)
        .finally(() => setPreparingDownload(false))
    } else {
      FileSaver.saveAs(downloadUrl, fileName)
    }
  }

  const Component = component

  const button = (
    <Component
      className={classes.root}
      disabled={disabled || preparingDownload} color={color} size={size}
      onClick={handleClick}
    >
      {children}
    </Component>
  )

  if (tooltip && !disabled && !preparingDownload) {
    return <Tooltip title={tooltip}>{button}</Tooltip>
  } else {
    return button
  }
})
Download.propTypes = {
  fileName: PropTypes.string,
  url: PropTypes.string,
  component: PropTypes.any,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]).isRequired,
  disabled: PropTypes.bool,
  tooltip: PropTypes.string,
  color: PropTypes.string,
  size: PropTypes.string
}

export default Download
