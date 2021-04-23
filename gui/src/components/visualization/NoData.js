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
import {
  Box,
  makeStyles
} from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'

/**
 * Component that is used as a replacement for data that is not available.
 * Notice that this is different from invalid data (this should display an
 * ErrorCard) or from a placeholder!
 */

// These styles do not depend on any props: they can be created once and are
// shared by each instance.
const useStaticStyles = makeStyles(theme => ({
  root: {
  },
  container: {
    backgroundColor: '#f3f3f3',
    width: '100%',
    height: '100%',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center'
  },
  placeholder: {
    position: 'absolute',
    top: theme.spacing(2),
    left: theme.spacing(2),
    right: theme.spacing(2),
    bottom: theme.spacing(2)
  },
  message: {
    color: theme.palette.primary.main,
    fontStyle: 'italic',
    opacity: 0.5,
    userSelect: 'none'
  }
}))
export default function NoData({aspectRatio, className, classes, 'data-testid': testID}) {
  // These styles change for each component individually
  const useStyles = makeStyles(theme => ({
    outerContainer: {
      height: 0,
      overflow: 'hidden',
      paddingBottom: `${100 / aspectRatio}%`,
      position: 'relative'
    }
  }))
  const staticStyles = useStaticStyles({classes: classes})
  const styles = useStyles()
  const content = <Box className={staticStyles.container} data-testid={testID}>
    <div className={staticStyles.message}>no data</div>
  </Box>
  return aspectRatio
    ? <div className={clsx(className, staticStyles.root)}>
      <div className={styles.outerContainer}>
        <div className={staticStyles.placeholder}>
          {content}
        </div>
      </div>
    </div>
    : content
}

NoData.propTypes = {
  aspectRatio: PropTypes.number,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}
