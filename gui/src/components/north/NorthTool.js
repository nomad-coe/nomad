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
import React, { useCallback, useContext, useEffect, useMemo, useState } from 'react'
import PropTypes from 'prop-types'
import { makeStyles } from '@material-ui/core/styles'
import {
  Box,
  Typography,
  Button
} from '@material-ui/core'
import AssessmentIcon from '@material-ui/icons/Assessment'
import Icon from '@material-ui/core/Icon'
import { useApi } from '../api.js'
import { useErrors } from '../errors'
import { northBase, enableNewHubApi } from '../../config'

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

export function getIconUrl(icon) {
  return (icon.startsWith('http://') || icon.startsWith('https://'))
    ? icon
    : `${process.env.PUBLIC_URL}/${icon}`
}

const LaunchButton = React.memo(function LaunchButton(props) {
  return <Button color="primary" variant="contained" size="small" {...props} />
})

export function useNorthToolHook(tool, uploadId, path) {
  const { name, path_prefix, with_path } = tool
  const { api } = useApi()
  const { raiseError } = useErrors()

  const getToolStatus = useCallback(() => {
    if (enableNewHubApi) {
      const extractMountState = response => {
        const uploadIds = response.uploads_ids
        if (uploadIds && typeof uploadIds === 'object') {
          return uploadId in uploadIds
        }
        return false
      }
      return api
        .get(`north/servers/${name}`)
        .then(response => ({ state: response.state, uploadid_is_mounted: extractMountState(response) }))
        .catch(error => {
          raiseError(error)
          return { state: "error", uploadid_is_mounted: null }
        })
    } else {
      return api
        .get(`north/${name}?upload_id=${uploadId}`)
        .then(response => ({ state: response.data.state, uploadid_is_mounted: response.upload_id_is_mounted }))
        .catch(error => {
          raiseError(error)
          return { state: "error", uploadid_is_mounted: null }
        })
    }
  }, [api, raiseError, name, uploadId])

  const launch = useCallback(() => {
    if (enableNewHubApi) {
      return api.post(`north/servers/${name}`)
        .then(() => {
          let toolUrl = `${northBase}/user-redirect/${name}`
          if (with_path && uploadId && path) {
            const params = new URLSearchParams()
            params.append('upload_id', uploadId)
            params.append('path', path)
            toolUrl = `${toolUrl}?${params.toString()}`
          }
          return { toolUrl, state: "running", uploadid_is_mounted: null }
        })
        .catch(error => {
          raiseError(error)
          return { toolUrl: null, state: "stopped", uploadid_is_mounted: null }
        })
    } else {
      return api.post(`north/${name}?upload_id=${uploadId}`)
        .then(response => {
          let toolUrl = `${northBase}/user/${response.username}/${response.tool}`

          if (with_path && response.upload_mount_dir && path) {
            if (path_prefix) {
              toolUrl = `${toolUrl}/${path_prefix}`
            }
            toolUrl = `${toolUrl}/${response.upload_mount_dir}/${path}`
          }

          return { toolUrl, state: response.data.state, uploadid_is_mounted: response.upload_id_is_mounted }
        })
        .catch(error => {
          raiseError(error)
          return { toolUrl: null, state: "stopped", uploadid_is_mounted: null }
        })
    }
  }, [api, raiseError, name, uploadId, path, path_prefix, with_path])

  const stop = useCallback(() => {
    const url = enableNewHubApi
      ? `north/servers/${name}`
      : `north/${name}?upload_id=${uploadId}`

    return api.delete(url)
      .then(() => "stopped")
      .catch(error => {
        raiseError(error)
        return "error"
      })
  }, [api, raiseError, name, uploadId])

  return { launch, stop, getToolStatus }
}

export const NorthToolButtons = React.memo(function NorthToolButton() {
  const { name, launch, stop, state } = useNorthTool()
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

const NorthTool = React.memo(function NorthTool({ tool, uploadId, path, children }) {
  const styles = useStyles()
  const [toolState, setToolState] = useState("stopped")
  const { launch, stop, getToolStatus } = useNorthToolHook(tool, uploadId, path)

  useEffect(() => {
    async function fetchStatus() {
      const { state: status } = await getToolStatus()
      setToolState(status)
    }
    fetchStatus()
  }, [getToolStatus])

  const handleLaunch = useCallback(async () => {
    setToolState("starting")
    const { toolUrl, state } = await launch()
    setToolState(state)
    if (toolUrl) window.open(toolUrl, tool.name)
  }, [launch, tool.name])

  const handleStop = useCallback(async () => {
    setToolState("stopping")
    const updatedState = await stop()
    setToolState(updatedState)
  }, [stop])

  const value = useMemo(() => ({
    state: toolState,
    launch: handleLaunch,
    stop: handleStop,
    ...tool
  }), [tool, toolState, handleLaunch, handleStop])

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
          {tool.icon ? (
            <Icon classes={{root: styles.iconRoot}}>
              <img className={styles.imageIcon} src={getIconUrl(tool.icon)} alt="icon"/>
            </Icon>
          ) : (
            <AssessmentIcon classes={{ root: styles.iconRoot }} />
          )}
          <Box flexGrow={1}>
            <Typography>
              <b>{tool.title || tool.name}{tool.version && <span> ({tool.version})</span>}</b>: {tool.short_description || tool.description}
            </Typography>
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
