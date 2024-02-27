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

import React, { useCallback, useState } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { Paper, Tooltip, IconButton } from '@material-ui/core'
import SearchIcon from '@material-ui/icons/Search'
import { useStyles } from '../search/SearchBar'
import { InputMetainfoControlled } from '../search/input/InputMetainfo'

/**
 * This component shows a search bar with autocomplete functionality.
 */
const ArchiveSearchBar = React.memo(({options, group, onChange, className}) => {
  const [value, setValue] = useState(null)
  const styles = useStyles()

  const handleSelect = useCallback((key) => {
    const obj = options[key]
    setValue(obj?.primary)
    onChange && onChange(obj.url)
  }, [options, onChange])

  const handleAccept = useCallback((key) => {
    const obj = options[key]
    if (obj && onChange) onChange(obj.url)
  }, [options, onChange])

  return <Paper className={clsx(className, styles.root)}>
    <InputMetainfoControlled
      value={value}
      onChange={setValue}
      onSelect={handleSelect}
      onAccept={handleAccept}
      options={options}
      group={group}
      TextFieldProps={{
        variant: "outlined",
        placeholder: "Type your keyword here",
        size: "medium"
      }}
      InputProps={{
        classes: {
          notchedOutline: styles.notchedOutline
        },
        startAdornment: <Tooltip title="Add filter">
          <IconButton className={styles.iconButton} onClick={() => handleAccept(value)} aria-label="search">
            <SearchIcon />
          </IconButton>
        </Tooltip>
      }}
    />
  </Paper>
})

ArchiveSearchBar.propTypes = {
  options: PropTypes.objectOf(PropTypes.shape({
    primary: PropTypes.string,
    secondary: PropTypes.string,
    group: PropTypes.string,
    key: PropTypes.string,
    url: PropTypes.string
  })),
  group: PropTypes.bool,
  onChange: PropTypes.func,
  className: PropTypes.string
}

export default ArchiveSearchBar
