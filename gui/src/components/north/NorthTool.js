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
import React, {useCallback, useEffect, useMemo, useState} from 'react'
import PropTypes from 'prop-types'
import { makeStyles } from '@material-ui/core/styles'
import {
  Box,
  Typography,
  Button
} from '@material-ui/core'
import AssessmentIcon from '@material-ui/icons/Assessment'
import Icon from '@material-ui/core/Icon'
import {useApi} from '../api.js'
import { useErrors } from '../errors'
import { northBase } from '../../config'

const launchButtonLabels = {
  'idle': 'Launch',
  'launching': 'Launching...',
  'running': 'Open',
  'stopping': 'Launch'
}

const stoppButtonLabels = {
  'stopping': 'Stopping...',
  'running': 'Stop'
}

const LaunchButton = React.memo(function LaunchButton(props) {
  return <Button color="primary" variant="contained" size="small" {...props} />
})

/**
 * Remote tool item. Can be shown within a list or as a standalone component.
 */
const useStyles = makeStyles(theme => ({
  root: {},
  imageIcon: {
    height: '100%'
  },
  iconRoot: {
    textAlign: 'center',
    height: '3rem',
    width: '3rem',
    marginRight: theme.spacing(2)
  },
  action: {
    marginTop: 0,
    marginRight: 0
  }
}))

const NorthTool = React.memo(({
  name,
  title,
  version,
  description,
  path_prefix,
  icon,
  disableDescription,
  disableActions,
  uploadId,
  path
}) => {
  const styles = useStyles()
  const {northApi, user} = useApi()
  const {raiseError} = useErrors()

  const [state, setState] = useState('idle')

  const toolUrl = useMemo(() => {
    let toolPath = ''
    if (path_prefix) {
      toolPath += `/${path_prefix}`
    }
    if (uploadId) {
      toolPath += `/uploads/${uploadId}`
    }
    if (path) {
      toolPath += `/${path}`
    }
    const toolUrl = `${northBase}/user/${user.preferred_username}/${name}${toolPath}`
    console.log('###', toolUrl)
    return toolUrl
  }, [user, name, path_prefix, uploadId, path])

  const getToolStatus = useCallback(() => {
    if (northApi === null) {
      return
    }
    return northApi.get(`servers/${name}/progress`)
      .then((response) => {
        const data = JSON.parse(response.data.substr(6))
        if (data.ready) {
          return 'running'
        } else {
          return 'launching'
        }
      })
      .catch(error => {
        if (error?.response?.status === 404 || error?.response?.status === 400) {
          return 'idle'
        } else {
          setState('idle')
          raiseError(error)
        }
      })
  }, [northApi, raiseError, name])

  useEffect(() => {
    getToolStatus().then((toolStatus) => {
      setState(toolStatus)
    })
  }, [setState, getToolStatus])

  const launch = useCallback(() => {
    // We get the current actual tools status and do not use the one used to display the status!
    getToolStatus().then((toolStatus) => {
      if (toolStatus === 'running') {
        window.open(toolUrl, name)
        setState(toolStatus)
      } else {
        setState('launching')
        northApi.post(`servers/${name}`)
          .then((response) => {
            window.open(toolUrl, name)
            setState('running')
          })
          .catch(errors => {
            raiseError(errors)
            setState('idle')
          })
      }
    })
  }, [setState, northApi, raiseError, name, toolUrl, getToolStatus])

  const stop = useCallback(() => {
    setState('stopping')
    northApi.delete(`servers/${name}`)
      .then((response) => {
        console.log(response)
        setState('idle')
      })
      .catch(raiseError)
  }, [northApi, raiseError, setState, name])

  return <Box marginY={1}>
    <Box display="flex" flexDirection="row" marginBottom={1}>
      {icon ? (
        <Icon classes={{root: styles.iconRoot}}>
          <img className={styles.imageIcon} src={`${process.env.PUBLIC_URL}/${icon}`} alt="icon"/>
        </Icon>
      ) : (
        <AssessmentIcon classes={{root: styles.iconRoot}}/>
      )}
      <Box flexGrow={1}>
        <Typography><b>{title || name}{version && <span> ({version})</span>}</b>: {description}</Typography>
      </Box>
    </Box>
    <Box display="flex" flexDirection="row">
      <LaunchButton fullWidth name={name} onClick={launch} disabled={state === 'stopping' || state === 'launching'}>
        {launchButtonLabels[state]}
      </LaunchButton>
      {(state === 'running' || state === 'stopping') && (
        <Box marginLeft={1}>
          <LaunchButton color="default" fullWidth onClick={stop} disabled={state === 'stopping'}>
            {stoppButtonLabels[state]}
          </LaunchButton>
        </Box>
      )}
    </Box>
  </Box>
})

NorthTool.propTypes = {
  name: PropTypes.string.isRequired, // The unique identitifier for this tool
  title: PropTypes.string.isRequired, // The displayed name of this tool
  version: PropTypes.string, // Version number if available
  description: PropTypes.string, // Version number if available
  path_prefix: PropTypes.string, // A path prefix. With path prefix selected upload or file will be added to path.
  icon: PropTypes.string, // Path to tool icon if available, relative to gui/public
  running: PropTypes.bool, // Whether the tool is running
  disableDescription: PropTypes.bool, // Whether to hide the description
  disableActions: PropTypes.bool, // Whether to hide actions for this tool
  uploadId: PropTypes.string, // Upload id
  path: PropTypes.string // path to a file within the upload
}

export default NorthTool
