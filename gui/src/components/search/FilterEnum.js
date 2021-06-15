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
import React, { useCallback, useEffect, useState } from 'react'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import {
  Grid,
  FormControlLabel,
  Checkbox
} from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { useApi } from '../apiV1'
import searchQuantities from '../../searchQuantities'
import FilterLabel from './FilterLabel'
import { useFilterState } from './FilterContext'

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    display: 'flex',
    alignItems: 'flex-start',
    justifyContent: 'center',
    flexDirection: 'column',
    boxSizing: 'border-box'
  },
  textField: {
    marginTop: theme.spacing(1)
  },
  input: {
    // padding: '16px 12px'
  }
}))
const FilterEnum = React.memo(({
  label,
  quantity,
  description,
  className,
  classes,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useStyles({classes: classes, theme: theme})
  const [options, setOptions] = useState()
  const api = useApi()

  // Determine the description and units
  const def = searchQuantities[quantity]
  const desc = description || def?.description || ''
  const name = label || def?.name
  const title = name

  // Attach the filter hook
  const [filter, setFilter] = useFilterState(quantity)

  // Fetch all possible values
  useEffect(() => {
    const search = {
      owner: 'visible',
      aggregations: {[quantity]: {terms: {quantity: quantity, size: 10}}},
      pagination: {page_size: 0},
      required: {include: []}
    }
    api.queryEntry(search)
      .then(data => {
        const opt = {}
        for (let value of data.aggregations[quantity].terms.data) {
          opt[value.value] = false
        }
        setOptions(opt)
      })
  }, [api, quantity])

  // React to changing filters
  useEffect(() => {
    if (filter) {
      setOptions(old => {
        const newOptions = {}
        for (let key of Object.keys(old)) {
          newOptions[key] = false
        }
        for (let value of filter) {
          newOptions[value] = true
        }
        return newOptions
      })
    }
  }, [filter])

  const handleChange = useCallback((event) => {
    // setOptions(old => {
    //   const newOptions = { ...old, [event.target.name]: event.target.checked }
    //   const checked = Object.entries(newOptions)
    //     .filter(([key, value]) => value)
    //     .map(([key, value]) => key)
    //   setFilter(checked)
    //   return newOptions
    // })

    const newOptions = { ...options, [event.target.name]: event.target.checked }
    const checked = Object.entries(newOptions)
      .filter(([key, value]) => value)
      .map(([key, value]) => key)
    setFilter(checked)
  }, [setFilter, options])

  const checkboxes = options && Object.entries(options).map(([key, value]) => {
    return <Grid item xs={12} key={key}>
      <FormControlLabel
        control={<Checkbox checked={value} onChange={handleChange} name={key}/>}
        label={key}
      />
    </Grid>
  })

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <FilterLabel label={title} description={desc}/>
    <Grid container spacing={0}>
      {checkboxes}
    </Grid>
  </div>
})

FilterEnum.propTypes = {
  label: PropTypes.string,
  quantity: PropTypes.string,
  description: PropTypes.string,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

export default FilterEnum
