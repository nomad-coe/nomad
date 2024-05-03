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
import {Quantity} from '../units/Quantity'
import {Unit} from '../units/Unit'
import {useUnitContext} from '../units/UnitContext'
import {UnitSelect} from './NumberEditQuantity'
import {getDisplayLabel} from "../../utils"
import {useRecoilValue} from "recoil"
import {configState} from "../archive/ArchiveBrowser"

export const SliderEditQuantity = React.memo((props) => {
  const {quantityDef, value, onChange, minValue, maxValue, ...sliderProps} = props
  const config = useRecoilValue(configState)
  const label = getDisplayLabel(quantityDef, true, config?.showMeta)

  const {units} = useUnitContext()
  const defaultUnit = useMemo(() => quantityDef.unit && new Unit(quantityDef.unit), [quantityDef])
  const dimension = defaultUnit && defaultUnit.dimension()
  const [unit, setUnit] = useState(units[dimension]?.definition || quantityDef.unit)
  const minValueConverted = useMemo(() => {
    return unit
      ? new Quantity(minValue, quantityDef.unit).to(unit).value()
      : minValue
  }, [minValue, quantityDef, unit])
  const maxValueConverted = useMemo(() => {
    return unit
      ? new Quantity(maxValue, quantityDef.unit).to(unit).value()
      : maxValue
  }, [maxValue, quantityDef, unit])
  const sliderValue = useMemo(() => {
    return unit
      ? new Quantity(value || 0, quantityDef.unit).to(unit).value()
      : value || 0
  }, [unit, quantityDef, value])

  const handleChangeValue = useCallback((event, value) => {
    const convertedValue = unit ? new Quantity(value, unit).to(quantityDef.unit).value() : value
    if (onChange) {
      onChange(convertedValue)
    }
  }, [unit, onChange, quantityDef])

  return <FormControl fullWidth>
    <FormLabel>{label}</FormLabel>
    <Box display="flex" alignItems="center">
      <Slider
        value={sliderValue}
        min={minValueConverted}
        max={maxValueConverted}
        onChange={handleChangeValue}
        valueLabelDisplay={(!unit ? 'on' : 'off')}
        {...sliderProps}
      />
      {unit && (
        <UnitSelect dimension={dimension} unit={unit} onChange={setUnit}/>
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
