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
 * wrapper around the MUI Skeleton component that allows using an aspect ratio
 * to determine the size.
 *
 * Override the 'placeholder' CSS-class to control the size of the placeholder
 * with respect to the root component.
 */

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    height: '100%'
  },
  placeholder: {
    position: 'absolute',
    top: theme.spacing(1),
    left: theme.spacing(2),
    right: theme.spacing(2),
    bottom: theme.spacing(1)
  },
  skeleton: {
    width: '100%',
    height: '100%'
  }
}))
const Placeholder = React.memo(({
  aspectRatio,
  className,
  classes,
  'data-testid':
  testID,
  ...other
}) => {
  // If aspect ratio is provided, use it to determine width and height
  const useStylesDynamic = makeStyles(theme => {
    return {
      containerOuter: aspectRatio
        ? {
          height: 0,
          overflow: 'hidden',
          paddingBottom: `${100 / aspectRatio}%`,
          position: 'relative'
        }
        : {
          width: '100%',
          height: '100%',
          position: 'relative'
        }
    }
  })
  const theme = useTheme()
  const styles = useStyles({classes: classes, theme: theme})
  const stylesDynamic = useStylesDynamic()
  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <div className={stylesDynamic.containerOuter}>
      <div className={styles.placeholder}>
        <Skeleton className={styles.skeleton} {...other} />
      </div>
    </div>
  </div>
})

Placeholder.propTypes = {
  aspectRatio: PropTypes.number,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

export default Placeholder
