
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
import React, { useMemo, useCallback } from 'react'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import { Tooltip, Typography, TextField, MenuItem } from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import searchQuantities from '../../searchQuantities'

const useStaticStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    alignItems: 'center',
    margin: `${theme.spacing(1)}px 0`
  },
  select: {
    margin: '0 0.5rem',
    flex: '0 0 12rem',
    minWidth: '12rem'
  },
  name: {
    marginLeft: theme.spacing(1),
    minWidth: '6rem'
  },
  textField: {
    marginTop: theme.spacing(1)
  }
}))
const SelectQuery = React.memo(({
  label,
  quantity,
  description,
  options,
  className,
  classes,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useStaticStyles({classes: classes, theme: theme})

  // Determine the description and options
  const [finalDescription, finalOptions] = useMemo(() => {
    const def = searchQuantities[quantity]
    const desc = description || def?.description || ''
    const opt = options || def?.type?.type_data || []
    return [desc, opt]
  }, [quantity, description, options])

  // Callback for handling input change
  const handleChange = useCallback((event) => {
  }, [])

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <Tooltip title={finalDescription}>
      <Typography className={styles.name} variant="body1">
        {`${label}:`}
      </Typography>
    </Tooltip>
    <TextField
      select
      defaultValue={finalOptions[0]}
      className={styles.select}
      margin='dense'
      size='small'
      variant='outlined'
      onChange={handleChange}
    >
      {finalOptions.map((option) => (
        <MenuItem key={option} value={option}>
          {option}
        </MenuItem>
      ))}
    </TextField>
  </div>
})

SelectQuery.propTypes = {
  label: PropTypes.string,
  quantity: PropTypes.string,
  description: PropTypes.string,
  options: PropTypes.arrayOf(PropTypes.string),
  className: PropTypes.string,
  classes: PropTypes.object,
  units: PropTypes.object,
  'data-testid': PropTypes.string
}

export default SelectQuery
