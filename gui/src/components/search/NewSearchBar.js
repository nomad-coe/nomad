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
import clsx from 'clsx'
import { makeStyles } from '@material-ui/core/styles'
import SearchIcon from '@material-ui/icons/Search'
import {
  MenuItem,
  Select,
  Paper,
  InputBase,
  Divider
} from '@material-ui/core'
import IconButton from '@material-ui/core/IconButton'

const useStyles = makeStyles(theme => ({
  root: {
    padding: '2px 4px',
    display: 'flex',
    alignItems: 'center',
    height: '2.8rem'
  },
  input: {
    marginLeft: theme.spacing(1),
    flex: 1
  },
  iconButton: {
    padding: 10
  },
  divider: {
    height: 'calc(100% - 6px)'
  },
  searchType: {
    paddingLeft: theme.spacing(1.5),
    paddingRight: theme.spacing(1.5)
  },
  select: {
    fontSize: '0.85rem'
  }
}))

/**
 * This searchbar component shows a searchbar with autocomplete functionality. The
 * searchbar also includes a status line about the current results. It uses the
 * search context to manipulate the current query and display results. It does its on
 * API calls to provide autocomplete suggestion options.
 */
const NewSearchBar = React.memo(({
  searchType,
  className,
  onSearchTypeChanged
}) => {
  const styles = useStyles()

  return <Paper className={clsx(className, styles.root)}>
    <div className={styles.searchType}>
      <Select
        value={searchType}
        onChange={onSearchTypeChanged}
        classes={{select: styles.select}}
      >
        <MenuItem value="nomad">NOMAD</MenuItem>
        <MenuItem value="optimade">OPTIMADE</MenuItem>
      </Select>
    </div>
    <Divider className={styles.divider} orientation="vertical" />
    <InputBase
      className={styles.input}
      placeholder={searchType === 'nomad' ? 'Search with quantity=value' : 'Search with Optimade filter language'}
      inputProps={{ 'aria-label': 'search' }}
    />
    <IconButton className={styles.iconButton} aria-label="search">
      <SearchIcon />
    </IconButton>
  </Paper>
})

NewSearchBar.propTypes = {
  searchType: PropTypes.string,
  className: PropTypes.string,
  onSearchTypeChanged: PropTypes.func
}

export default NewSearchBar
