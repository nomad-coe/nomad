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
import React, {useCallback, useState} from 'react'
import {Box, FormControlLabel, Collapse, TextField, MenuItem} from '@material-ui/core'
import PropTypes from 'prop-types'
import ArrowDownIcon from '@material-ui/icons/ArrowDropDown'
import ArrowRightIcon from '@material-ui/icons/ArrowRight'
import {NumberFieldWithUnit, useNumberEditQuantityStyles} from './NumberEditQuantity'
import {Unit, useUnits} from '../../units'
import {conversionMap, unitMap} from '../../unitsData'
import {StringField, WithHelp, getFieldProps} from './StringEditQuantity'

export const ListEditQuantity = React.memo((props) => {
  const {quantityDef, section, component, componentProps, defaultValues, onChange, actions} = props
  const [values, setValues] = useState(defaultValues)
  const [open, setOpen] = useState(false)
  let Component = component

  const handleChange = useCallback((newValue, index) => {
    let newValues = [...values]
    newValues[index] = newValue
    setValues(newValues)
    if (onChange) {
      onChange(newValues, section, quantityDef)
    }
  }, [onChange, quantityDef, section, values])

  if (!defaultValues) return ''

  return <Box display='block'>
    <Box display='flex'>
      <WithHelp helpDescription={quantityDef.description} {...getFieldProps(quantityDef)}>
        <Box sx={{ flexDirection: 'column' }} onClick={() => setOpen(!open)} height={'26px'}>
          <FormControlLabel
            control={open ? <ArrowDownIcon/> : <ArrowRightIcon/>}
            label={getFieldProps(quantityDef).label}
          />
        </Box>
      </WithHelp>
      {open && actions}
    </Box>
    <Collapse in={open}>
      <div>
        {defaultValues.map((value, index) => {
          return <Box key={index} marginTop={1}>
            <Component
              defaultValue={value !== undefined ? value : ''} onChange={newValue => handleChange(newValue, index)} {...componentProps} label={undefined}
              inputProps={{style: { padding: '14px' }}}
            />
          </Box>
        })}
      </div>
    </Collapse>
  </Box>
})
ListEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  component: PropTypes.any.isRequired,
  componentProps: PropTypes.object.isRequired,
  defaultValues: PropTypes.arrayOf(PropTypes.any),
  onChange: PropTypes.func.isRequired,
  actions: PropTypes.element
}

export const ListStringEditQuantity = React.memo((props) => {
  const {quantityDef, section, onChange, ...otherProps} = props
  const shape = quantityDef.type?.shape
  const defaultValue = (quantityDef.default !== undefined ? quantityDef.default : Array.apply(null, Array(shape[0])).map(() => ''))
  let values = section[quantityDef.name] || defaultValue

  return <ListEditQuantity
    quantityDef={quantityDef}
    section={section}
    component={StringField}
    componentProps={otherProps}
    defaultValues={values}
    onChange={onChange}/>
})
ListStringEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

export const ListNumberEditQuantity = React.memo((props) => {
  const classes = useNumberEditQuantityStyles()
  const {quantityDef, section, onChange, ...otherProps} = props
  const systemUnits = useUnits()
  const shape = quantityDef.type?.shape
  const defaultValue = (quantityDef.default !== undefined ? quantityDef.default : Array.apply(null, Array(shape[0])).map(() => ''))
  const dimension = quantityDef.unit && unitMap[quantityDef.unit].dimension
  const units = quantityDef.unit && conversionMap[dimension].units
  const isUnit = quantityDef.unit && ['float64', 'float32', 'float'].includes(quantityDef.type?.type_data)
  const [unit, setUnit] = useState(systemUnits[dimension] || quantityDef.unit)
  let values = section[quantityDef.name] || defaultValue

  const handleChangeUnit = useCallback((newUnit) => {
    setUnit(newUnit)
  }, [])

  const componentProps = {
    defaultUnit: quantityDef.unit,
    dataType: quantityDef.type?.type_data,
    unit: unit,
    ...otherProps
  }

  return <ListEditQuantity
    quantityDef={quantityDef}
    section={section}
    component={NumberFieldWithUnit}
    componentProps={componentProps}
    defaultValues={values}
    onChange={onChange}
    actions={isUnit && <TextField
      className={classes.unitSelect} variant='filled' size='small' select
      label="unit" value={unit}
      onChange={(event) => handleChangeUnit(event.target.value)}
    >
      {units.map(unit => <MenuItem key={unit} value={unit}>{(new Unit(unit)).label()}</MenuItem>)}
    </TextField>}/>
})
ListNumberEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}
