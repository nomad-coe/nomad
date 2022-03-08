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
import React, {useCallback, useMemo, useState} from 'react'
import {
  Box,
  FormControl,
  FormLabel, Slider
} from '@material-ui/core'
import PropTypes from 'prop-types'
import {convertUnit, useUnits} from '../../units'
import {unitMap} from '../../unitsData'
import {UnitSelect} from './NumberEditQuantity'
import {getFieldProps} from './StringEditQuantity'

export const SliderEditQuantity = React.memo((props) => {
  const {quantityDef, value, onChange, minValue, maxValue, ...sliderProps} = props
  const {label} = getFieldProps(quantityDef)

  const systemUnits = useUnits()
  const hasUnit = quantityDef.unit
  const dimension = hasUnit && unitMap[quantityDef.unit].dimension
  const [unit, setUnit] = useState(systemUnits[dimension] || quantityDef.unit)

  const sliderValue = useMemo(() => {
    if (!hasUnit) {
      return value || 0
    }
    return convertUnit(value || 0, quantityDef.unit, unit)
  }, [hasUnit, unit, quantityDef, value])

  const handleChangeValue = useCallback((event, value) => {
    const convertedValue = hasUnit ? convertUnit(value, unit, quantityDef.unit) : value
    if (onChange) {
      onChange(convertedValue)
    }
  }, [hasUnit, unit, onChange, quantityDef])

  return <FormControl fullWidth>
    <FormLabel>{label}</FormLabel>
    <Box display="flex" alignItems="center">
      <Slider
        value={sliderValue}
        min={convertUnit(minValue, quantityDef.unit, unit)}
        max={convertUnit(maxValue, quantityDef.unit, unit)}
        onChange={handleChangeValue}
        valueLabelDisplay={(!hasUnit ? 'on' : 'off')}
        {...sliderProps}
      />
      {hasUnit && (
        <UnitSelect defaultUnit={quantityDef.unit} unit={unit} onChange={setUnit}/>
      )}
    </Box>
  </FormControl>
})
SliderEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.number,
  minValue: PropTypes.number.isRequired,
  maxValue: PropTypes.number.isRequired,
  onChange: PropTypes.func
}
