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
import React, {useEffect, useState} from 'react'
import PropTypes from 'prop-types'
import { makeStyles } from '@material-ui/core/styles'
import {
  CardHeader,
  CardContent,
  Typography,
  Button
} from '@material-ui/core'
import Grid from '@material-ui/core/Grid';
import AssessmentIcon from '@material-ui/icons/Assessment'
import Icon from '@material-ui/core/Icon'
import NORTHLaunchButton from './NORTHLaunchButton'
import { useApi } from '../api'

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
  },
  stopButton: {
    color: 'white',
    backgroundColor: '#FF2D00',
    '&:hover': {
      backgroundColor: '#CC2400',
    },
    width: '100%'
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

  const path = path_prefix && uploadId ? `${path_prefix}/uploads/${uploadId}` : null

  const {api, user} = useApi()

  const [launching, setLaunching] = useState(false)
  const [running, setRunning] = useState(false)

  useEffect(() => {
    api.axios.get(`http://localhost:9000/fairdi/nomad/latest/north/hub/api/users/${user.preferred_username}/servers/${name}/progress`, { headers:{'Authorization': 'token 51dc0e836aee4d51a10681b970943e3d'} })
    .then((response) => {
      console.log(JSON.parse(response.data.substr(response.data.indexOf('{'))).ready)
      setRunning(JSON.parse(response.data.substr(response.data.indexOf('{'))).ready)
    })
    .catch((error) => {
      console.log(error)
    })
  }, [user.preferred_username])

  const tryLaunch = async () => {
    api.axios.get(`http://localhost:9000/fairdi/nomad/latest/north/hub/api/users/${user.preferred_username}/servers/${name}/progress`, { headers:{'Authorization': 'token 51dc0e836aee4d51a10681b970943e3d'} })
    .then((response) => {
      console.log(JSON.parse(response.data.substr(response.data.indexOf('{'))).url)
      window.open('http://localhost:9000' + JSON.parse(response.data.substr(response.data.indexOf('{'))).url, '_blank').focus();
      setLaunching(false)
      setRunning(true)
    })
    .catch((error) => {
      console.log(error)
      if (error.response)
        console.log(error.response.status)
        setLaunching(true)
        api.axios.post(`http://localhost:9000/fairdi/nomad/latest/north/hub/api/users/${user.preferred_username}/servers/${name}`, {}, { headers:{'Authorization': 'token 51dc0e836aee4d51a10681b970943e3d'} })
        .then((response) => {
          tryLaunch()
        })
    })
  }

  const stop = () => {
    api.axios.delete(`http://localhost:9000/fairdi/nomad/latest/north/hub/api/users/${user.preferred_username}/servers/${name}`, { headers:{'Authorization': 'token 51dc0e836aee4d51a10681b970943e3d'} })
    .then((response) => {
      console.log(response)
      setRunning(false)
    })
  }

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
            <NORTHLaunchButton
              name={name}
              path={path}
              onClick={tryLaunch}
            >{launching ? 'Launching...' : (running ? "Open" : "Launch")}</NORTHLaunchButton>
          </Grid>
          <Grid item xs={12}>
            {running && <Button variant="contained" className={styles.stopButton} onClick={stop}>Stop</Button>}
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
