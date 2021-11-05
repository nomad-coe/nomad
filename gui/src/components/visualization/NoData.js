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
import { Typography, makeStyles } from '@material-ui/core'
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
    height: '100%'
  },
  relative: {
    width: '100%',
    height: '100%',
    position: 'relative'
  },
  placeholder: {
    position: 'absolute',
    top: theme.spacing(0),
    left: theme.spacing(1),
    right: theme.spacing(1),
    bottom: theme.spacing(0),
    backgroundColor: theme.palette.grey[100],
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    borderRadius: theme.spacing(0.5)
  },
  message: {
    color: theme.palette.primary.main,
    fontStyle: 'italic',
    opacity: 0.5,
    userSelect: 'none'
  }
}))
export default function NoData({className, classes, 'data-testid': testID}) {
  const styles = useStyles({classes: classes})

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <div className={styles.relative}>
      <div className={styles.placeholder}>
        <Typography className={styles.message}>no data</Typography>
      </div>
    </div>
  </div>
}

NoData.propTypes = {
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}
