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
import React, { useMemo } from 'react'
import { isArray } from 'lodash'
import { Typography, makeStyles, useTheme } from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'

/**
 * Component that is used as a replacement for data that is not available.
 * Notice that this is different from invalid data (this should display an
 * ErrorCard) or from a placeholder!
 */

// These styles do not depend on any props: they can be created once and are
// shared by each instance.
const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    height: '100%',
    position: 'relative'
  },
  placeholder: {
    width: '100%',
    height: '100%',
    backgroundColor: theme.palette.grey[100],
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center'
  },
  message: {
    color: theme.palette.primary.main,
    fontStyle: 'italic',
    opacity: 0.5,
    userSelect: 'none'
  }
}))
export default function NoData({
  className,
  margin,
  'data-testid': testID
}) {
  const styles = useStyles()
  const theme = useTheme()
  const margins = useMemo(() => {
    let marginTop, marginRight, marginBottom, marginLeft
    if (isArray(margin)) {
      marginTop = margin[0]
      marginRight = margin[1]
      marginBottom = margin[2]
      marginLeft = margin[3]
    } else {
      marginTop = margin
      marginRight = margin
      marginBottom = margin
      marginLeft = margin
    }
    return {
      position: 'absolute',
      top: theme.spacing(marginTop),
      left: theme.spacing(marginLeft),
      right: theme.spacing(marginRight),
      bottom: theme.spacing(marginBottom)
    }
  }, [theme, margin])

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <div style={margins}>
      <div className={styles.placeholder}>
        <Typography className={styles.message}>no data</Typography>
      </div>
    </div>
  </div>
}

NoData.propTypes = {
  className: PropTypes.string,
  /**
   * Margin with respect to the parent in theme spacing units.
   */
  margin: PropTypes.oneOfType([PropTypes.number, PropTypes.arrayOf(PropTypes.number)]),
  'data-testid': PropTypes.string
}

NoData.defaultProps = {
  margin: [0, 1, 0, 1]
}
