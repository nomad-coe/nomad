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
import React, {useCallback, useEffect, useState} from 'react'
import {
  TextField,
  Box,
  MenuItem,
  FormControl,
  FormLabel, Slider
} from '@material-ui/core'
import PropTypes from 'prop-types'
import {convertUnit, Unit, useUnits} from '../../units'
import {conversionMap, unitMap} from '../../unitsData'
import {useNumberEditQuantityStyles} from './NumberEditQuantity'
import {getFieldProps} from './StringEditQuantity'

export const SliderEditQuantity = React.memo((props) => {
  const classes = useNumberEditQuantityStyles()
  const {quantityDef, section, onChange, ...otherProps} = props
  const {minValue, maxValue, ...sliderProps} = otherProps
  const defaultValue = (quantityDef.default !== undefined ? quantityDef.default : undefined)
  const [value, setValue] = useState(0)
  const [convertedValue, setConvertedValue] = useState(0)
  const dimension = quantityDef.unit && unitMap[quantityDef.unit].dimension
  const units = quantityDef.unit && conversionMap[dimension].units
  const systemUnits = useUnits()
  const isUnit = quantityDef.unit && ['float64', 'float32', 'float'].includes(quantityDef.type?.type_data)
  const [unit, setUnit] = useState(systemUnits[dimension] || quantityDef.unit)

  useEffect(() => {
    let newValue = section[quantityDef.name] !== undefined ? section[quantityDef.name] : (defaultValue || minValue)
    setValue(newValue)
    setConvertedValue(`${(isUnit ? (!isNaN(Number(newValue)) || newValue === '' ? convertUnit(Number(newValue), quantityDef.unit, unit) : '') : newValue)}`)
  }, [defaultValue, isUnit, minValue, quantityDef, section, unit])

  const handleChangeUnit = useCallback((newUnit) => {
    setUnit(newUnit)
    setConvertedValue(`${(isUnit ? (!isNaN(Number(value)) || value === '' ? convertUnit(Number(value), quantityDef.unit, newUnit) : '') : value)}`)
  }, [isUnit, quantityDef, value])

  const handleChangeValue = useCallback((event, newValue) => {
    if (typeof newValue !== 'number') return
    setConvertedValue(`${newValue}`)
    if (onChange) {
      onChange((isUnit ? (newValue === '' ? newValue : (!isNaN(Number(newValue)) ? convertUnit(Number(newValue), unit, quantityDef.unit) : '')) : newValue), section, quantityDef)
    }
    setValue((isUnit ? (!isNaN(Number(newValue)) || newValue === '' ? convertUnit(Number(newValue), unit, quantityDef.unit) : '') : newValue))
    setConvertedValue(`${Number(newValue)}`)
  }, [isUnit, unit, onChange, quantityDef, section])

  return <FormControl fullWidth>
    <FormLabel>{getFieldProps(quantityDef).label}</FormLabel>
    <Box display='flex'>
      <Slider
        value={Number(convertedValue)}
        min={convertUnit(Number(minValue), quantityDef.unit, unit)}
        max={convertUnit(Number(maxValue), quantityDef.unit, unit)}
        onChange={handleChangeValue}
        valueLabelDisplay={(!isUnit ? 'on' : 'off')}
        {...sliderProps}/>
      {isUnit && <TextField
        className={classes.unitSelect} variant='filled' size='small' select
        label="unit" value={unit}
        onChange={(event) => handleChangeUnit(event.target.value)}
      >
        {units.map(unit => <MenuItem key={unit} value={unit}>{(new Unit(unit)).label()}</MenuItem>)}
      </TextField>}
    </Box>
  </FormControl>
})
SliderEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}
