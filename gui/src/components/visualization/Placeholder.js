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
import { makeStyles, useTheme } from '@material-ui/core/styles'
import { Skeleton } from '@material-ui/lab'
import PropTypes from 'prop-types'
import clsx from 'clsx'

/**
 * Component that is used as a placeholder while loading data. Fairly simple
 * wrapper around the MUI Skeleton component.
 */

// These styles do not depend on any props: they can be created once and are
// shared by each instance.
const useStaticStyles = makeStyles(theme => ({
  root: {
  },
  placeholder: {
    position: 'absolute',
    top: theme.spacing(2),
    left: theme.spacing(2),
    right: theme.spacing(2),
    bottom: theme.spacing(2)
  },
  skeleton: {
    width: '100%',
    height: '100%'
  }
}))
export default function Placeholder(props) {
  // If aspect ratio is provided, use it to determine width and height
  const {aspectRatio, className, classes, ...other} = props
  const useStyles = makeStyles(theme => {
    const style = {}
    if (aspectRatio) {
      style.containerOuter = {
        height: 0,
        overflow: 'hidden',
        paddingBottom: `${100 / aspectRatio}%`,
        position: 'relative'
      }
    }
    return style
  })
  const theme = useTheme()
  const styles = useStyles()
  const staticStyles = useStaticStyles({classes: classes, theme: theme})
  if (aspectRatio) {
    return <div className={clsx(className, styles.root)}>
      <div className={styles.containerOuter}>
        <div className={staticStyles.placeholder}>
          <Skeleton variant="rect" className={staticStyles.skeleton} {...other}/>
        </div>
      </div>
    </div>
  }
  return <div className={clsx(className, staticStyles.root)}>
    <div className={styles.containerInner}>
      <Skeleton {...other} className={staticStyles.skeleton}></Skeleton>
    </div>
  </div>
}

Placeholder.propTypes = {
  aspectRatio: PropTypes.number,
  className: PropTypes.string,
  classes: PropTypes.object
}
