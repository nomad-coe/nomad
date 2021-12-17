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
import {
  Grid,
  Divider,
  makeStyles
} from '@material-ui/core'

/**
 * For displaying a grid of input properties.
 */
const inputGridSpacing = 2
export function InputGrid({children}) {
  return <Grid container spacing={inputGridSpacing} style={{marginTop: 0}}>
    {children}
  </Grid>
}

InputGrid.propTypes = {
  children: PropTypes.any
}

/**
 * For displaying an individual input filter, typically within a InputGrid.
 */
const useInputGridItemStyles = makeStyles(theme => ({
  root: {
    marginBottom: theme.spacing(inputGridSpacing / 2)
  },
  divider: {
    marginLeft: theme.spacing(-inputGridSpacing),
    marginRight: theme.spacing(-inputGridSpacing),
    marginTop: theme.spacing(-inputGridSpacing / 2),
    marginBottom: theme.spacing(inputGridSpacing / 2),
    backgroundColor: theme.palette.grey[300]
  }
}))
export function InputGridItem({classes, children, ...other}) {
  const styles = useInputGridItemStyles({classes: classes})
  return <Grid item {...other} className={styles.root}>
    <Divider className={styles.divider}/>
    {children}
  </Grid>
}

InputGridItem.propTypes = {
  classes: PropTypes.object,
  children: PropTypes.any
}

InputGridItem.defaultProps = {
}
