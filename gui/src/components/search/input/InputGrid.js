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
import React, { useContext } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import {
  Grid,
  Divider,
  makeStyles
} from '@material-ui/core'
import { inputSectionContext } from './InputNestedObject'

/**
 * For displaying a grid of input properties.
 *
 * Notice that we have to use a workaround for the MUI grid item in order to
 * prevent horizontal overflow. A padding is applied to the parent element as
 * described in https://v4.mui.com/components/grid/#negative-margin
 */

/**
 * For displaying an individual input filter, typically within a InputGrid.
 */
export const inputItemPaddingHorizontal = 0.75
export const inputItemPaddingVertical = 1
const useStyles = makeStyles(theme => ({
  root: {},
  padding: {
    paddingLeft: theme.spacing(inputItemPaddingHorizontal),
    paddingRight: theme.spacing(inputItemPaddingHorizontal)
  }
}))

export function InputGrid({className, disablePadding, children}) {
  const styles = useStyles()
  return <div className={clsx(className, styles.root, !disablePadding && styles.padding)}>
    <Grid container spacing={0}>
      {children}
    </Grid>
  </div>
}

InputGrid.propTypes = {
  children: PropTypes.any,
  disablePadding: PropTypes.bool,
  className: PropTypes.string
}

/**
 * For displaying an individual input filter, typically within a InputGrid.
 */
const useInputGridItemStyles = makeStyles(theme => ({
  root: {
  },
  padding: {
    padding: theme.spacing(inputItemPaddingVertical, inputItemPaddingHorizontal, inputItemPaddingVertical, inputItemPaddingHorizontal)
  },
  divider: {
    backgroundColor: theme.palette.grey[300]
  },
  paddingEnabledDivider: {
    marginTop: theme.spacing(-inputItemPaddingVertical),
    marginRight: theme.spacing(-2 * inputItemPaddingHorizontal),
    marginLeft: theme.spacing(-2 * inputItemPaddingHorizontal)
  },
  paddingDisabled: {
    marginRight: theme.spacing(-inputItemPaddingHorizontal),
    marginLeft: theme.spacing(-inputItemPaddingHorizontal)
  },
  section: {
    marginTop: theme.spacing(-2)
  }
}))
export function InputGridItem({classes, children, disablePadding, ...other}) {
  const sectionContext = useContext(inputSectionContext)
  const section = sectionContext?.section
  const styles = useInputGridItemStyles({classes: classes})
  return <Grid item {...other} className={clsx(styles.root, !disablePadding && styles.padding)}>
    {section
      ? <div className={styles.section}></div>
      : <Divider className={clsx(styles.divider, disablePadding ? styles.paddingDisabled : styles.paddingEnabledDivider)}/>
    }
    <div className={clsx(disablePadding && styles.paddingDisabled)}>
      {children}
    </div>
  </Grid>
}

InputGridItem.propTypes = {
  classes: PropTypes.object,
  children: PropTypes.any,
  disablePadding: PropTypes.any
}

InputGridItem.defaultProps = {
  xs: 12
}
