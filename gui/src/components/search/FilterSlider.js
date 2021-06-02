
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
import React, { useState, useMemo, useCallback } from 'react'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import {
  Slider,
  TextField,
  FormHelperText
} from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import FilterLabel from './FilterLabel'
import { Quantity, Unit } from '../../units'
import searchQuantities from '../../searchQuantities'
import { useSetFilter } from './FilterContext'

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
    width: '6rem'
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
const FilterText = React.memo(({
  label,
  quantity,
  description,
  step,
  min,
  max,
  start,
  end,
  className,
  classes,
  units,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useStyles({classes: classes, theme: theme})
  const setFilter = useSetFilter(quantity)
  // const {stats, subscribe, unsubscribe} = useStatistics(quantity)
  // console.log(stats)

  // // Subscribe to stats on mount, unsubscribe on unmounting
  // useEffect(() => {
  //   subscribe({quantity: quantity})
  //   return () => unsubscribe()
  // // eslint-disable-next-line react-hooks/exhaustive-deps
  // }, [quantity])

  // TODO: get the minimum and maximum from aggregations
  const {trueMin, trueMax} = useMemo(() => {
    return {
      trueMin: 1,
      trueMax: 10
    }
  }, [])
  const [trueStart, setTrueStart] = useState(start || trueMin)
  const [trueEnd, setTrueEnd] = useState(end || trueMax)
  const [range, setRange] = useState({gte: trueStart, lte: trueEnd})
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

  // Handle start change: whenever a valid number is set, the value is committed
  // to the hook.
  const handleStartChange = useCallback((event) => {
    const value = event.target?.value
    setTrueStart(value)
    const number = Number.parseFloat(value)
    if (!isNaN(number)) {
      const outOfRange = number < trueMin
      if (outOfRange) {
        setError(`Minimum value cannot be below ${trueMin}.`)
        return
      }
      setError()
      setRange(old => {
        const newRange = {...old, gte: number}
        sendFilter(newRange)
        return newRange
      })
    }
  }, [sendFilter, trueMin])

  // Handle end change: whenever a valid number is set, the value is committed
  // to the hook.
  const handleEndChange = useCallback((event) => {
    const value = event.target?.value
    setTrueEnd(value)
    const number = Number.parseFloat(value)
    if (!isNaN(number)) {
      const outOfRange = number > trueMax
      if (outOfRange) {
        setError(`Maximum value cannot be below ${trueMax}.`)
        return
      }
      setError()
      setRange(old => {
        const newRange = {...old, lte: number}
        sendFilter(newRange)
        return newRange
      })
    }
  }, [sendFilter, trueMax])

  // Handle range commit: Set the filter when mouse is released on a slider
  const handleRangeCommit = useCallback((event, value) => {
    sendFilter({gte: value[0], lte: value[1]})
  }, [sendFilter])

  // Handle range change: only change the rendered values, send to the filter
  // hook only after mouseup
  const handleRangeChange = useCallback((event, value) => {
    setTrueStart(value[0])
    setTrueEnd(value[1])
    setRange({gte: value[0], lte: value[1]})
    setError()
  }, [])

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <FilterLabel label={title} description={desc}/>
    <div className={styles.inputRow}>
      <TextField
        variant="outlined"
        label="min"
        className={styles.textField}
        value={trueStart}
        margin="normal"
        onChange={handleStartChange}
        InputProps={{classes: {input: styles.input}}}
      />
      <div className={styles.spacer}>
        <Slider
          min={trueMin}
          max={trueMax}
          step={1}
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
        value={trueEnd}
        margin="normal"
        onChange={handleEndChange}
        InputProps={{classes: {input: styles.input}}}
      />
    </div>
    {error && <FormHelperText error>
      {error}
    </FormHelperText>}
  </div>
})

FilterText.propTypes = {
  label: PropTypes.string,
  quantity: PropTypes.string,
  description: PropTypes.string,
  min: PropTypes.number,
  max: PropTypes.number,
  start: PropTypes.number,
  end: PropTypes.number,
  step: PropTypes.number,
  className: PropTypes.string,
  classes: PropTypes.object,
  units: PropTypes.object,
  'data-testid': PropTypes.string
}

export default FilterText
