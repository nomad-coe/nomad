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
import React, {useCallback, useContext, useEffect, useMemo, useState} from 'react'
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

const northToolContext = React.createContext()

export function useNorthTool() {
  return useContext(northToolContext)
}

const launchButtonLabels = {
  'stopped': 'Launch',
  'starting': 'Launching...',
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

export const NorthToolButtons = React.memo(function NorthToolButton() {
  const {name, launch, stop, state} = useNorthTool()
  return (
    <Box display="flex" flexDirection="row">
      <LaunchButton fullWidth name={name} onClick={launch} disabled={state === 'stopping' || state === 'starting' || !state}>
        {launchButtonLabels[state] || 'not available'}
      </LaunchButton>
      {(state === 'running' || state === 'stopping') && (
        <Box marginLeft={1}>
          <LaunchButton color="default" fullWidth onClick={stop} disabled={state === 'stopping'}>
            {stoppButtonLabels[state]}
          </LaunchButton>
        </Box>
      )}
    </Box>
  )
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

const NorthTool = React.memo(function NorthTool({tool, uploadId, path, children}) {
  const {name, title, version, description, short_description, icon, path_prefix, with_path} = tool
  const styles = useStyles()
  const {api} = useApi()
  const {raiseError} = useErrors()

  const [state, setState] = useState('stopped')

  const getToolStatus = useCallback(() => {
    return api.get(`north/${name}`)
      .then(response => {
        return response.data.state
      }).catch(raiseError)
  }, [api, raiseError, name])

  useEffect(() => {
    const toolStatus = getToolStatus()
    if (toolStatus) {
      toolStatus.then((toolStatus) => {
        setState(toolStatus)
      })
    }
  }, [setState, getToolStatus])

  const launch = useCallback(() => {
    // We get the current actual tools status and do not use the one used to display the status!
    setState('starting')
    api.post(`north/${name}?upload_id=${uploadId}`)
      .then((response) => {
        console.log(response)
        let toolUrl = `${northBase}/user/${response.username}/${response.tool}`
        if (with_path && response.upload_mount_dir && path) {
          if (path_prefix) {
            toolUrl = `${toolUrl}/${path_prefix}`
          }
          toolUrl = `${toolUrl}/${response.upload_mount_dir}/${path}`
        }
        console.log(toolUrl)
        window.open(toolUrl, name)
        setState(response.data.state)
      })
      .catch(errors => {
        raiseError(errors)
        setState('stopped')
      })
  }, [setState, api, raiseError, name, with_path, uploadId, path, path_prefix])

  const stop = useCallback(() => {
    setState('stopping')
    api.delete(`north/${name}`)
      .then((response) => {
        console.log(response)
        setState('stopped')
      })
      .catch(raiseError)
  }, [api, raiseError, setState, name])

  const value = useMemo(() => ({
    state: state,
    launch: launch,
    stop: stop,
    ...tool
  }), [tool, state, launch, stop])

  if (children) {
    return (
      <northToolContext.Provider value={value}>
        {children}
      </northToolContext.Provider>
    )
  }

  return (
    <northToolContext.Provider value={value}>
      <Box marginY={1}>
        <Box display="flex" flexDirection="row" marginBottom={1}>
          {icon ? (
            <Icon classes={{root: styles.iconRoot}}>
              <img className={styles.imageIcon} src={`${process.env.PUBLIC_URL}/${icon}`} alt="icon"/>
            </Icon>
          ) : (
            <AssessmentIcon classes={{root: styles.iconRoot}}/>
          )}
          <Box flexGrow={1}>
            <Typography><b>{title || name}{version && <span> ({version})</span>}</b>: {short_description || description}</Typography>
          </Box>
        </Box>
        <NorthToolButtons />
      </Box>
    </northToolContext.Provider>
  )
})

NorthTool.propTypes = {
  tool: PropTypes.object.isRequired,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]),
  uploadId: PropTypes.string, // the upload id with the path in it
  path: PropTypes.string // path to a file within the upload
}

export default NorthTool
