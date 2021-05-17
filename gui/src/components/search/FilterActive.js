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
import { makeStyles } from '@material-ui/core/styles'
import {
  Paper,
  IconButton,
  Tooltip,
} from '@material-ui/core'
import ClearIcon from '@material-ui/icons/Clear'
import { searchContext } from './SearchContext'

/**
 * Displays the active filters.
 */
const useStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    alignItems: 'center',
    padding: '2px 4px'
  },
  filterIcon: {
    paddingLeft: theme.spacing(0.5),
    paddingRight: theme.spacing(0.5)
  },
  spacer: {
    flexGrow: 1
  },
  divider: {
    height: 'calc(100% - 6px)'
  }
}))
const FiltersActive = React.memo(({
  className
}) => {
  const styles = useStyles()
  const {query} = useContext(searchContext)

  return <Paper className={clsx(className, styles.root)}>
    <div className={styles.spacer}></div>
    <Tooltip
      title="Clear filters"
    >
      <IconButton className={styles.iconButton} aria-label="search">
        <ClearIcon />
      </IconButton>
    </Tooltip>
  </Paper>
})
FiltersActive.propTypes = {
  className: PropTypes.string
}

export default FiltersActive
