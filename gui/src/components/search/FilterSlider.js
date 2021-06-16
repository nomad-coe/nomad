
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
import React, { useState, useMemo, useCallback, useEffect } from 'react'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import {
  Slider,
  TextField,
  FormHelperText
} from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import FilterLabel from './FilterLabel'
import { Quantity, Unit, toUnitSystem } from '../../units'
import searchQuantities from '../../searchQuantities'
import { useFilterState, useAgg } from './FilterContext'
import { formatNumber } from '../../utils'

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
    marginTop: theme.spacing(1),
    flexGrow: 1,
    width: '7rem'
  },
  spacer: {
    flex: '2 1 100%',
    paddingLeft: '18px',
    paddingRight: '18px'
  },
  inputRow: {
    width: '100%',
    display: 'flex',
    justifyContent: 'space-between',
    alignItems: 'center'
  },
  thumb: {
    '&:active': {
      boxShadow: '0px 0px 0px 12px rgb(0, 141, 195, 16%)'
    },
    '&:focusVisible': {
      boxShadow: '0px 0px 0px 6px rgb(0, 141, 195, 16%)'
    }
  },
  input: {
    padding: '16px 12px'
  }
}))
const FilterSlider = React.memo(({
  label,
  quantity,
  description,
  step,
  visible,
  className,
  classes,
  units,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useStyles({classes: classes, theme: theme})
  const [filter, setFilter] = useFilterState(quantity)
  const agg = useAgg(quantity, 'min_max', true, visible)
  const [minGlobal, maxGlobal] = agg || [undefined, undefined]
  const [minText, setMinText] = useState()
  const [maxText, setMaxText] = useState()
  const [minLocal, setMinLocal] = useState()
  const [maxLocal, setMaxLocal] = useState()
  const [range, setRange] = useState()
  const [error, setError] = useState()

  // Determine the description and units
  const def = searchQuantities[quantity]
  const desc = description || def?.description || ''
  const name = label || def?.name
  const unitSI = def?.unit
  const unit = useMemo(() => {
    return unitSI && new Unit(unitSI, units)
  }, [unitSI, units])
  const unitLabel = unit && unit.label()
  const title = unitLabel ? `${name} (${unitLabel})` : name
  const stepConverted = step instanceof Quantity ? step.toSystem(units) : step
  let minConverted = (minGlobal !== undefined && unitSI) ? toUnitSystem(minGlobal, unitSI, units) : minGlobal
  let maxConverted = (maxGlobal !== undefined && unitSI) ? toUnitSystem(maxGlobal, unitSI, units) : maxGlobal

  // If not manual range has been specified, the range is automatically adjusted
  // according to global min/max of the field
  useEffect(() => {
    if (maxConverted !== undefined && minConverted !== undefined) {
      setMaxLocal(maxConverted)
      setMinLocal(minConverted)
      if (filter === undefined) {
        setMaxText(formatNumber(maxConverted))
        setMinText(formatNumber(minConverted))
        setRange({gte: minConverted, lte: maxConverted})
      }
    }
  }, [maxConverted, minConverted, filter])

  // Function for converting search values and sending them to the search
  // context.
  const sendFilter = useCallback(range => {
    if (unit) {
      range = {
        lte: new Quantity(range.lte, unit),
        gte: new Quantity(range.gte, unit)
      }
    }
    setFilter(range)
  }, [unit, setFilter])

  const handleStartChange = useCallback((event) => {
    setError()
    const value = event.target?.value
    setMinText(value)
  }, [])

  const handleEndChange = useCallback((event) => {
    setError()
    const value = event.target?.value
    setMaxText(value)
  }, [])

  const handleStartSubmit = useCallback(() => {
    const value = minText
    const number = Number.parseFloat(value)
    if (!isNaN(number)) {
      const outOfRange = number < minConverted
      if (outOfRange) {
        setError(`Minimum value cannot be below ${minConverted}.`)
        return
      }
      setRange(old => {
        const newRange = {...old, gte: number}
        sendFilter(newRange)
        return newRange
      })
    } else {
      setError(`Invalid value.`)
    }
  }, [minText, minConverted, sendFilter])

  const handleEndSubmit = useCallback(() => {
    const value = maxText
    const number = Number.parseFloat(value)
    if (!isNaN(number)) {
      const outOfRange = number > maxConverted
      if (outOfRange) {
        setError(`Maximum value cannot be above ${maxConverted}.`)
        return
      }
      setRange(old => {
        const newRange = {...old, gte: number}
        sendFilter(newRange)
        return newRange
      })
    } else {
      setError(`Invalid value.`)
    }
  }, [maxText, maxConverted, sendFilter])

  // Handle range commit: Set the filter when mouse is released on a slider
  const handleRangeCommit = useCallback((event, value) => {
    sendFilter({gte: value[0], lte: value[1]})
  }, [sendFilter])

  // Handle range change: only change the rendered values, send to the filter
  // hook only after mouseup
  const handleRangeChange = useCallback((event, value) => {
    setMinText(value[0])
    setMaxText(value[1])
    setRange({gte: value[0], lte: value[1]})
    setError()
  }, [])

  // Until an initial min/max range is available, we do not display the
  // component
  if (range === undefined) {
    return
  }

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <FilterLabel label={title} description={desc}/>
    <div className={styles.inputRow}>
      <TextField
        variant="outlined"
        label="min"
        className={styles.textField}
        value={minText}
        margin="normal"
        onChange={handleStartChange}
        onBlur={handleStartSubmit}
        onKeyDown={(event) => { if (event.key === 'Enter') { handleStartSubmit() } }}
        InputProps={{classes: {input: styles.input}}}
      />
      <div className={styles.spacer}>
        <Slider
          min={minLocal}
          max={maxLocal}
          step={stepConverted}
          value={[range.gte, range.lte]}
          onChange={handleRangeChange}
          onChangeCommitted={handleRangeCommit}
          valueLabelDisplay="off"
          classes={{thumb: styles.thumb, active: styles.active}}
        />
      </div>
      <TextField
        variant="outlined"
        label="max"
        className={styles.textField}
        value={maxText}
        margin="normal"
        onChange={handleEndChange}
        onBlur={handleEndSubmit}
        onKeyDown={(event) => { if (event.key === 'Enter') { handleEndSubmit() } }}
        InputProps={{classes: {input: styles.input}}}
      />
    </div>
    {error && <FormHelperText error>
      {error}
    </FormHelperText>}
  </div>
})

FilterSlider.propTypes = {
  label: PropTypes.string,
  quantity: PropTypes.string.isRequired,
  description: PropTypes.string,
  step: PropTypes.oneOfType([PropTypes.number, PropTypes.object]).isRequired,
  visible: PropTypes.bool,
  className: PropTypes.string,
  classes: PropTypes.object,
  units: PropTypes.object,
  'data-testid': PropTypes.string
}

export default FilterSlider
