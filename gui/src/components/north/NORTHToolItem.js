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
import { makeStyles } from '@material-ui/core/styles'
import {
  CardHeader,
  CardContent,
  Typography
} from '@material-ui/core'
import AssessmentIcon from '@material-ui/icons/Assessment'
import Icon from '@material-ui/core/Icon'
import NORTHLaunchButton from './NORTHLaunchButton'

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
  running,
  disableDescription,
  disableActions,
  uploadId
}) => {
  const styles = useStyles()

  const path = path_prefix && uploadId ? `${path_prefix}/uploads/${uploadId}` : null

  return <div>
    <CardHeader
      avatar={icon
        ? <Icon classes={{root: styles.iconRoot}}>
          <img className={styles.imageIcon} src={`../${icon}`} alt="icon"/>
        </Icon>
        : <AssessmentIcon classes={{root: styles.iconRoot}}/>}
      title={title}
      subheader={version}
      action={!disableActions && <NORTHLaunchButton
        name={name}
        path={path}
      >Launch</NORTHLaunchButton>}
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
