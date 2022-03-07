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
import { makeStyles, useTheme } from '@material-ui/core/styles'
import { Skeleton } from '@material-ui/lab'
import PropTypes from 'prop-types'
import clsx from 'clsx'

/**
 * Component that is used as a placeholder while loading data. Implemented as a
 * simple wrapper around the MUI Skeleton component.
 *
 * Override the 'placeholder' CSS-class to control the size of the placeholder
 * with respect to the root component.
 */

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    height: '100%',
    position: 'relative'
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
  className,
  classes,
  'data-testid':
  testID,
  margin,
  ...other
}) => {
  const styles = useStyles({classes: classes})
  const theme = useTheme()
  const margins = useMemo(() => {
    return {
      position: 'absolute',
      top: theme.spacing(margin),
      left: theme.spacing(margin),
      right: theme.spacing(margin),
      bottom: theme.spacing(margin)
    }
  }, [margin, theme])

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <div style={margins}>
      <Skeleton variant='rect' className={styles.skeleton} {...other} />
    </div>
  </div>
})

Placeholder.propTypes = {
  className: PropTypes.string,
  classes: PropTypes.object,
  /**
   * Placeholder margin with respect to the parent in theme spacing units.
   */
  margin: PropTypes.oneOfType([PropTypes.number, PropTypes.arrayOf(PropTypes.number)]),
  'data-testid': PropTypes.string
}

Placeholder.defaultProps = {
  margin: 0
}

export default Placeholder
