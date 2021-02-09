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
export default function NoData({aspectRatio, className, classes}) {
  const useStyles = makeStyles(theme => {
    return {
      root: {
      },
      outerContainer: {
        height: 0,
        overflow: 'hidden',
        paddingBottom: `${100 / aspectRatio}%`,
        position: 'relative'
      },
      innerContainer: {
        position: 'absolute',
        top: 0,
        left: 0,
        width: '100%',
        height: '100%'
      },
      box: {
        width: '100%',
        height: '100%',
        boxSizing: 'border-box',
        padding: theme.spacing(2)
      },
      background: {
        backgroundColor: '#f3f3f3',
        width: '100%',
        height: '100%',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center'
      },
      message: {
        color: theme.palette.primary.main,
        opacity: 0.5,
        userSelect: 'none'
      }
    }
  })
  const styles = useStyles(classes)
  const content = <Box className={styles.box}>
    <Box className={styles.background}>
      <div className={styles.message}>NO DATA</div>
    </Box>
  </Box>
  return aspectRatio
    ? <div className={clsx(className, styles.root)}>
      <div className={styles.outerContainer}>
        <div className={styles.innerContainer}>
          {content}
        </div>
      </div>
    </div>
    : content
}

NoData.propTypes = {
  aspectRatio: PropTypes.number,
  className: PropTypes.string,
  classes: PropTypes.string
}
