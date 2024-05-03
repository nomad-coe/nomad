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
import PropTypes from 'prop-types'
import { number, bool, reach } from 'yup'
import { cloneDeep } from 'lodash'
import {
  Checkbox,
  FormControlLabel
} from '@material-ui/core'
import { InputJMESPath } from '../input/InputMetainfo'
import { schemaWidget, schemaAxis, schemaMarkers } from './Widget'
import { WidgetEditDialog, WidgetEditGroup, WidgetEditOption, WidgetEditSelect } from './WidgetEdit'
import { useSearchContext } from '../SearchContext'
import { autorangeDescription } from './WidgetHistogram'
import { DType, setDeep, parseJMESPath } from '../../../utils'
import { InputTextField } from '../input/InputText'
import UnitInput from '../../units/UnitInput'

// Predefined in order to not break memoization
const dtypesNumeric = new Set([DType.Int, DType.Float, DType.Timestamp])
const dtypesColor = new Set([DType.String, DType.Enum, DType.Float, DType.Int])
const nPointsOptions = {
  100: 100,
  1000: 1000,
  10000: 10000
}
function isEmptyString(value) {
  return value === undefined || value === null || !value?.trim?.()?.length
}
/**
 * A dialog that is used to configure a scatter plot widget.
 */
export const WidgetScatterPlotEdit = React.memo(({widget}) => {
    const { filterData, useSetWidget } = useSearchContext()
    const [settings, setSettings] = useState(cloneDeep(widget))
    const [errors, setErrors] = useState({})
    const [dimensions, setDimensions] = useState({})
    const setWidget = useSetWidget(widget.id)

    const handleError = useCallback((key, value) => {
      setErrors(old => ({...old, [key]: value}))
    }, [setErrors])

    const handleErrorQuantity = useCallback((key, value) => {
      handleError(key, value)
      setDimensions((old) => ({...old, [key]: null}))
    }, [handleError])

    const handleChange = useCallback((key, value) => {
      setSettings(old => {
        const newValue = {...old}
        setDeep(newValue, key, value)
        return newValue
      })
    }, [setSettings])

    const handleClose = useCallback(() => {
      setWidget(old => ({...old, editing: false}))
    }, [setWidget])

    const handleAccept = useCallback((key, value) => {
      try {
        reach(schemaWidgetScatterPlot, key).validateSync(value)
      } catch (e) {
        handleError(key, e.message)
        return
      }
      setErrors(old => ({...old, [key]: undefined}))
      handleChange(key, value)
    }, [handleError, handleChange])

    const handleAcceptQuantity = useCallback((key, value) => {
      handleAccept(key, value)
      const { quantity } = parseJMESPath(value)
      const dimension = filterData[quantity]?.dimension
      setDimensions((old) => ({...old, [key]: dimension}))
    }, [handleAccept, filterData])

    // Upon accepting the entire form, we perform final validation that also
    // takes into account cross-field incompatibilities
    const handleEditAccept = useCallback(() => {
      // Check for independent errors from components
      const independentErrors = Object.values(errors).some(x => !!x)
      if (independentErrors) return

      // Check for missing values: TODO: This kind of check should be replaced
      // by a dedicated form context. Inputs could automatically register into
      // this form to enable imperative validation etc.
      const xEmpty = isEmptyString(settings?.x?.quantity)
      if (xEmpty) {
        handleErrorQuantity('x.quantity', 'Please specify a value')
      }
      const yEmpty = isEmptyString(settings?.y?.quantity)
      if (yEmpty) {
        handleErrorQuantity('y.quantity', 'Please specify a value')
      }

      if (!independentErrors && !xEmpty && !yEmpty) {
        setWidget(old => ({...old, ...{...settings, editing: false, visible: true}}))
      }
    }, [settings, setWidget, errors, handleErrorQuantity])

    return <WidgetEditDialog
        id={widget.id}
        open={widget.editing}
        visible={widget.visible}
        title="Edit scatter plot widget"
        onClose={handleClose}
        onAccept={handleEditAccept}
      >
      <WidgetEditGroup title="x axis">
        <WidgetEditOption>
          <InputJMESPath
            label="quantity"
            value={settings.x?.quantity}
            onChange={(value) => handleChange('x.quantity', value)}
            onSelect={(value) => handleAcceptQuantity('x.quantity', value)}
            onAccept={(value) => handleAcceptQuantity('x.quantity', value)}
            error={errors['x.quantity']}
            onError={(value) => handleErrorQuantity('x.quantity', value)}
            dtypes={dtypesNumeric}
            dtypesRepeatable={dtypesNumeric}
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <InputTextField
            label="title"
            fullWidth
            value={settings.x?.title}
            onChange={(event) => handleChange('x.title', event.target.value)}
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <UnitInput
            label='unit'
            value={settings.x?.unit}
            onChange={(value) => handleChange('x.unit', value)}
            onSelect={(value) => handleAccept('x.unit', value)}
            onAccept={(value) => handleAccept('x.unit', value)}
            error={errors['x.unit']}
            onError={(value) => handleError('x.unit', value)}
            dimension={dimensions['x.quantity'] || null}
            optional
            disableGroup
          />
        </WidgetEditOption>
      </WidgetEditGroup>
      <WidgetEditGroup title="y axis">
        <WidgetEditOption>
          <InputJMESPath
            label="quantity"
            value={settings.y?.quantity}
            onChange={(value) => handleChange('y.quantity', value)}
            onSelect={(value) => handleAcceptQuantity('y.quantity', value)}
            onAccept={(value) => handleAcceptQuantity('y.quantity', value)}
            error={errors['y.quantity']}
            onError={(value) => handleErrorQuantity('y.quantity', value)}
            dtypes={dtypesNumeric}
            dtypesRepeatable={dtypesNumeric}
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <InputTextField
            label="title"
            fullWidth
            value={settings.y?.title}
            onChange={(event) => handleChange('y.title', event.target.value)}
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <UnitInput
            label='unit'
            value={settings.y?.unit}
            onChange={(value) => handleChange('y.unit', value)}
            onSelect={(value) => handleAccept('y.unit', value)}
            onAccept={(value) => handleAccept('y.unit', value)}
            error={errors['y.unit']}
            onError={(value) => handleError('y.unit', value)}
            dimension={dimensions['y.quantity'] || null}
            optional
            disableGroup
          />
        </WidgetEditOption>
      </WidgetEditGroup>
      <WidgetEditGroup title="marker color">
        <WidgetEditOption>
          <InputJMESPath
            label="quantity"
            value={settings?.markers?.color?.quantity}
            onChange={(value) => handleChange('markers.color.quantity', value)}
            onSelect={(value) => handleAcceptQuantity('markers.color.quantity', value)}
            onAccept={(value) => handleAcceptQuantity('markers.color.quantity', value)}
            error={errors['markers.color.quantity']}
            onError={(value) => handleErrorQuantity('markers.color.quantity', value)}
            dtypes={dtypesColor}
            dtypesRepeatable={dtypesColor}
            optional
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <InputTextField
            label="title"
            fullWidth
            value={settings.markers?.color?.title}
            onChange={(event) => handleChange('markers.color.title', event.target.value)}
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <UnitInput
            label='unit'
            value={settings.markers?.color?.unit}
            onChange={(value) => handleChange('markers.color.unit', value)}
            onSelect={(value) => handleAccept('markers.color.unit', value)}
            onAccept={(value) => handleAccept('markers.color.unit', value)}
            error={errors['markers.color.unit']}
            onError={(value) => handleError('markers.color.unit', value)}
            dimension={dimensions['markers.color.quantity'] || null}
            optional
            disableGroup
          />
        </WidgetEditOption>
      </WidgetEditGroup>
      <WidgetEditGroup title="general">
        <WidgetEditOption>
          <InputTextField
            label="title"
            fullWidth
            value={settings?.title}
            onChange={(event) => handleChange('title', event.target.value)}
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <WidgetEditSelect
            label="Maximum number of entries to load"
            options={nPointsOptions}
            value={settings.size}
            onChange={(event) => { handleChange('size', event.target.value) }}
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <FormControlLabel
            control={<Checkbox checked={settings.autorange} onChange={(event, value) => handleChange('autorange', value)}/>}
            label={autorangeDescription}
          />
        </WidgetEditOption>
      </WidgetEditGroup>
    </WidgetEditDialog>
})

WidgetScatterPlotEdit.propTypes = {
  widget: PropTypes.object
}

export const schemaWidgetScatterPlot = schemaWidget.shape({
  x: schemaAxis.required('Quantity for the x axis is required.'),
  y: schemaAxis.required('Quantity for the y axis is required.'),
  markers: schemaMarkers,
  size: number().integer().required('Size is required.'),
  autorange: bool()
})
