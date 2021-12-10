
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
import React, { useState, useMemo, useCallback, useEffect, useRef } from 'react'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import { Slider, FormHelperText } from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { isNil } from 'lodash'
import InputHeader from './InputHeader'
import InputTooltip from './InputTooltip'
import { InputTextField } from './InputText'
import { Quantity, Unit, toUnitSystem, toSI } from '../../../units'
import { formatNumber } from '../../../utils'
import searchQuantities from '../../../searchQuantities'
import { useSearchContext } from '../SearchContext'

function format(value) {
  return formatNumber(value, 'float', 6, true)
}

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
    marginTop: 0,
    marginBotton: 0,
    flexGrow: 1,
    width: '10rem'
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
  }
}))
const InputSlider = React.memo(({
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
  const {filterData, useAgg, useFilterState, useFilterLocked} = useSearchContext()
  const styles = useStyles({classes: classes, theme: theme})
  const endChanged = useRef(false)
  const startChanged = useRef(false)
  const [filter, setFilter] = useFilterState(quantity)
  const locked = useFilterLocked(quantity)
  const agg = useAgg(quantity, visible)
  const [minGlobalSI, maxGlobalSI] = agg?.data || [undefined, undefined]
  const [minText, setMinText] = useState('')
  const [maxText, setMaxText] = useState('')
  const [minLocal, setMinLocal] = useState(0)
  const [maxLocal, setMaxLocal] = useState(0)
  const [range, setRange] = useState({gte: undefined, lte: undefined})
  const [error, setError] = useState()

  // Determine the description and units
  const def = searchQuantities[quantity]
  const desc = description || def?.description || ''
  const name = label || def?.name
  const unitSI = filterData[quantity].unit || def?.unit || 'dimensionless'
  const unit = useMemo(() => {
    return unitSI && new Unit(unitSI, units)
  }, [unitSI, units])
  const unitLabel = unit && unit.label()
  const title = unitLabel ? `${name} (${unitLabel})` : name
  const stepSI = step instanceof Quantity ? step.toSI() : step
  const unavailable = (minGlobalSI === null || maxGlobalSI === null || range === undefined)
  const disabled = locked || unavailable

  // The slider minimum and maximum are set according to global min/max of the
  // field.
  useEffect(() => {
    if (!isNil(minGlobalSI) && !isNil(maxGlobalSI)) {
      setMinLocal(minGlobalSI)
      setMaxLocal(maxGlobalSI)
    }
  }, [minGlobalSI, maxGlobalSI])

  // When units change or the slider is used to set a value, update the min/max
  // text
  useEffect(() => {
    setMinText(isNil(range.gte) ? '' : format(toUnitSystem(range.gte, unitSI, units)))
    setMaxText(isNil(range.lte) ? '' : format(toUnitSystem(range.lte, unitSI, units)))
  }, [range, units, unitSI])

  // If no filter has been specified by the user, the range is automatically
  // adjusted according to global min/max of the field. If filter is set, the
  // slider value is set according to it.
  useEffect(() => {
    let lte
    let gte
    let min
    let max
    if (!isNil(minGlobalSI) && !isNil(maxGlobalSI)) {
      // When no filter is set, use the whole available range
      if (isNil(filter)) {
        gte = minGlobalSI
        lte = maxGlobalSI
        min = minGlobalSI
        max = maxGlobalSI
      // A single specific value is given
      } else if (filter instanceof Quantity) {
        gte = filter.toSI()
        lte = filter.toSI()
        min = Math.min(minGlobalSI, gte)
        max = Math.max(maxGlobalSI, lte)
      // A range is given
      } else {
        gte = filter.gte instanceof Quantity ? filter.gte.toSI() : filter.gte
        lte = filter.lte instanceof Quantity ? filter.lte.toSI() : filter.lte
        if (isNil(gte)) {
          min = Math.min(minGlobalSI, lte)
          gte = Math.min(minGlobalSI, lte)
        } else {
          min = Math.min(minGlobalSI, gte)
        }
        if (isNil(lte)) {
          max = Math.max(maxGlobalSI, gte)
          lte = Math.max(maxGlobalSI, gte)
        } else {
          max = Math.max(maxGlobalSI, lte)
        }
      }
      setMinLocal(min)
      setMaxLocal(max)
      setRange({gte: gte, lte: lte})
      setMinText(format(toUnitSystem(gte, unitSI, units)))
      setMaxText(format(toUnitSystem(lte, unitSI, units)))
    }
  }, [minGlobalSI, maxGlobalSI, filter, unitSI, units])

  // Function for converting search values and sending them to the search
  // context.
  const sendFilter = useCallback(range => {
    if (unitSI) {
      range = {
        lte: new Quantity(range.lte, unitSI),
        gte: new Quantity(range.gte, unitSI)
      }
    }
    setFilter(range)
  }, [unitSI, setFilter])

  // Used to simultaneously update the range in the slider and the actual filter
  // value.
  const updateRange = useCallback((range) => {
    setRange(range)
    sendFilter(range)
  }, [setRange, sendFilter])

  // Handles changes in the text input fields
  const handleChange = useCallback((ref, setter) => {
    return (event) => {
      ref.current = true
      setError()
      const value = event.target?.value
      setter(value)
    }
  }, [])

  // Called when min values are submitted through the text field.
  const handleMinSubmit = useCallback((a) => {
    if (!startChanged.current) {
      return
    }
    const value = minText
    const number = Number.parseFloat(value)
    if (!isNaN(number)) {
      const numberSI = toSI(number, unit)
      const outOfRange = numberSI < minGlobalSI
      if (outOfRange) {
        const min = toUnitSystem(minGlobalSI, unitSI, units)
        setError(`Minimum value cannot be below ${min}.`)
        return
      }
      updateRange({...range, gte: numberSI})
    } else {
      setError(`Invalid minimum value.`)
    }
    startChanged.current = false
  }, [range, minText, unitSI, unit, units, minGlobalSI, updateRange])

  // Called when max values are submitted through the text field.
  const handleMaxSubmit = useCallback(() => {
    if (!endChanged.current) {
      return
    }
    const value = maxText
    const number = Number.parseFloat(value)
    if (!isNaN(number)) {
      const numberSI = toSI(number, unit)
      const outOfRange = numberSI > maxGlobalSI
      if (outOfRange) {
        const max = toUnitSystem(maxGlobalSI, unitSI, units)
        setError(`Maximum value cannot be above ${max}.`)
        return
      }
      updateRange({...range, lte: numberSI})
    } else {
      setError(`Invalid maximum value.`)
    }
    endChanged.current = false
  }, [range, maxText, unitSI, unit, units, maxGlobalSI, updateRange])

  // Handle range commit: Set the filter when mouse is released on a slider
  const handleRangeCommit = useCallback((event, value) => {
    sendFilter({gte: value[0], lte: value[1]})
  }, [sendFilter])

  // Handle range change: only change the rendered values, send to the filter
  // hook only after mouseup
  const handleRangeChange = useCallback((event, value) => {
    setRange({gte: value[0], lte: value[1]})
    setError()
  }, [])

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <InputHeader
      quantity={quantity}
      label={title}
      description={desc}
      disableStatistics
      disableScale
    />
    <InputTooltip locked={locked} unavailable={unavailable}>
      <div>
        <div className={styles.inputRow}>
          <InputTextField
            disabled={disabled}
            label="min"
            className={styles.textField}
            value={minText}
            margin="normal"
            onChange={handleChange(startChanged, setMinText)}
            onBlur={handleMinSubmit}
            onKeyDown={(event) => { if (event.key === 'Enter') { handleMinSubmit() } }}
          />
          <div className={styles.spacer}>
            <Slider
              disabled={disabled}
              color="secondary"
              min={minLocal}
              max={maxLocal}
              step={stepSI}
              value={[range.gte, range.lte]}
              onChange={handleRangeChange}
              onChangeCommitted={handleRangeCommit}
              valueLabelDisplay="off"
              classes={{thumb: styles.thumb, active: styles.active}}
            />
          </div>
          <InputTextField
            disabled={disabled}
            label="max"
            className={styles.textField}
            value={maxText}
            margin="normal"
            onChange={handleChange(endChanged, setMaxText)}
            onBlur={handleMaxSubmit}
            onKeyDown={(event) => { if (event.key === 'Enter') { handleMaxSubmit() } }}
          />
        </div>
        {error && <FormHelperText error>
          {error}
        </FormHelperText>}
      </div>
    </InputTooltip>
  </div>
})

InputSlider.propTypes = {
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

export default InputSlider
