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
import { makeStyles } from '@material-ui/core/styles'
import { Skeleton } from '@material-ui/lab'
import PropTypes from 'prop-types'
import clsx from 'clsx'

/**
 * Component that is used as a placeholder while loading data. Fairly simple
 * wrapper around the MUI Skeleton component.
 */
export default function Placeholder(props) {
  // If aspect ratio is provided, use it to determine width and height
  const {aspectRatio, className, classes, ...other} = props
  const useStyles = makeStyles(props => {
    if (aspectRatio) {
      return {
        root: {
        },
        skeletonContainer: {
          height: 0,
          overflow: 'hidden',
          paddingBottom: `${100 / aspectRatio}%`,
          position: 'relative'
        },
        skeleton: {
          position: 'absolute',
          top: 0,
          left: 0,
          width: '100%',
          height: '100%'
        }
      }
    }
  })
  const styles = useStyles(classes)
  if (aspectRatio) {
    return <div className={clsx(className, styles.root)}>
      <div className={styles.skeletonContainer}>
        <Skeleton variant="rect" className={styles.skeleton} {...other}/>
      </div>
    </div>
  }
  return <Skeleton {...other}></Skeleton>
}

Placeholder.propTypes = {
  aspectRatio: PropTypes.number,
  className: PropTypes.string,
  classes: PropTypes.string
}
