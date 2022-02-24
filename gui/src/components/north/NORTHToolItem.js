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
  CardHeader,
  CardContent,
  Typography,
  Button
} from '@material-ui/core'
import Grid from '@material-ui/core/Grid'
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
  return <Button color="primary" style={{width: '8rem'}} variant="contained" {...props} />
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

const NORTHToolItem = React.memo(({
  name,
  title,
  version,
  description,
  path_prefix,
  icon,
  disableDescription,
  disableActions,
  uploadId
}) => {
  const styles = useStyles()
  const {northApi, user} = useApi()
  const {raiseError} = useErrors()

  const [state, setState] = useState('idle')

  const toolUrl = useMemo(() => {
    const path = path_prefix && uploadId ? `/${path_prefix}/uploads/${uploadId}` : ''
    return `${northBase}/user/${user.preferred_username}/${name}${path}`
  }, [user, name, path_prefix, uploadId])

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

  return <div>
    <CardHeader
      avatar={icon
        ? <Icon classes={{root: styles.iconRoot}}>
          <img className={styles.imageIcon} src={`../${icon}`} alt="icon"/>
        </Icon>
        : <AssessmentIcon classes={{root: styles.iconRoot}}/>}
      title={title}
      subheader={version}
      action={!disableActions &&
        <Grid container direction="column" spacing={1}>
          <Grid item xs={12}>
            <LaunchButton name={name} onClick={launch} disabled={state === 'stopping' || state === 'launching'}>
              {launchButtonLabels[state]}
            </LaunchButton>
          </Grid>
          <Grid item xs={12}>
            {(state === 'running' || state === 'stopping') && <LaunchButton color="error" onClick={stop} disabled={state === 'stopping'}>
              {stoppButtonLabels[state]}
            </LaunchButton>}
          </Grid>
        </Grid>
      }
      classes={{action: styles.action}}
    />
    {(!disableDescription && description) && <CardContent>
      <Typography variant="body2" component="p">{description}</Typography>
    </CardContent>}
  </div>
})

NORTHToolItem.propTypes = {
  name: PropTypes.string.isRequired, // The unique identitifier for this tool
  title: PropTypes.string.isRequired, // The displayed name of this tool
  version: PropTypes.string, // Version number if available
  description: PropTypes.string, // Version number if available
  path_prefix: PropTypes.string, // A path prefix. With path prefix selected upload or file will be added to path.
  icon: PropTypes.string, // Path to tool icon if available, relative to gui/public
  running: PropTypes.bool, // Whether the tool is running
  disableDescription: PropTypes.bool, // Whether to hide the description
  disableActions: PropTypes.bool, // Whether to hide actions for this tool
  uploadId: PropTypes.string // Upload id
}

export default NORTHToolItem
