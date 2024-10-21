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
  FormControlLabel,
  MenuItem,
  TextField
} from '@material-ui/core'
import { InputJMESPath } from '../input/InputMetainfo'
import { schemaWidget, schemaAxis, schemaMarkers } from './Widget'
import { WidgetEditDialog, WidgetEditGroup, WidgetEditOption, WidgetEditSelect } from './WidgetEdit'
import { useSearchContext } from '../SearchContext'
import { autorangeDescription } from './WidgetHistogram'
import { DType, setDeep, parseJMESPath, isEmptyString } from '../../../utils'
import { InputTextField } from '../input/InputText'
import { scalesLimited } from '../../plotting/common'
import UnitInput from '../../units/UnitInput'

// Predefined in order to not break memoization
const dtypesNumeric = new Set([DType.Int, DType.Float, DType.Timestamp])
const dtypesColor = new Set([DType.String, DType.Enum, DType.Float, DType.Int])
const nPointsOptions = {
  100: 100,
  1000: 1000,
  10000: 10000
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

      // Check for missing values. This check is required because there is no
      // value set when a new widget is created, and pressing the done button
      // without filling a value should raise an error.
      const xEmpty = isEmptyString(settings?.x?.search_quantity)
      if (xEmpty) {
        handleErrorQuantity('x.search_quantity', 'Please specify a value.')
      }
      const yEmpty = isEmptyString(settings?.y?.search_quantity)
      if (yEmpty) {
        handleErrorQuantity('y.search_quantity', 'Please specify a value.')
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
            label="Search quantity"
            value={settings.x?.search_quantity}
            onChange={(value) => handleChange('x.search_quantity', value)}
            onSelect={(value) => handleAcceptQuantity('x.search_quantity', value)}
            onAccept={(value) => handleAcceptQuantity('x.search_quantity', value)}
            error={errors['x.search_quantity']}
            onError={(value) => handleErrorQuantity('x.search_quantity', value)}
            dtypes={dtypesNumeric}
            dtypesRepeatable={dtypesNumeric}
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <InputTextField
            label="Title"
            fullWidth
            value={settings.x?.title}
            onChange={(event) => handleChange('x.title', event.target.value)}
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <UnitInput
            label='Unit'
            value={settings.x?.unit}
            onChange={(value) => handleChange('x.unit', value)}
            onSelect={(value) => handleAccept('x.unit', value)}
            onAccept={(value) => handleAccept('x.unit', value)}
            error={errors['x.unit']}
            onError={(value) => handleError('x.unit', value)}
            dimension={dimensions['x.search_quantity'] || null}
            optional
            disableGroup
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <TextField
            select
            fullWidth
            label="Scale"
            variant="filled"
            value={settings.x?.scale}
            onChange={(event) => { handleChange('x.scale', event.target.value) }}
          >
            {Object.keys(scalesLimited).map((key) =>
              <MenuItem value={key} key={key}>{key}</MenuItem>
            )}
          </TextField>
        </WidgetEditOption>
      </WidgetEditGroup>
      <WidgetEditGroup title="y axis">
        <WidgetEditOption>
          <InputJMESPath
            label="Search quantity"
            value={settings.y?.search_quantity}
            onChange={(value) => handleChange('y.search_quantity', value)}
            onSelect={(value) => handleAcceptQuantity('y.search_quantity', value)}
            onAccept={(value) => handleAcceptQuantity('y.search_quantity', value)}
            error={errors['y.search_quantity']}
            onError={(value) => handleErrorQuantity('y.search_quantity', value)}
            dtypes={dtypesNumeric}
            dtypesRepeatable={dtypesNumeric}
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <InputTextField
            label="Title"
            fullWidth
            value={settings.y?.title}
            onChange={(event) => handleChange('y.title', event.target.value)}
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <UnitInput
            label='Unit'
            value={settings.y?.unit}
            onChange={(value) => handleChange('y.unit', value)}
            onSelect={(value) => handleAccept('y.unit', value)}
            onAccept={(value) => handleAccept('y.unit', value)}
            error={errors['y.unit']}
            onError={(value) => handleError('y.unit', value)}
            dimension={dimensions['y.search_quantity'] || null}
            optional
            disableGroup
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <TextField
            select
            fullWidth
            label="Scale"
            variant="filled"
            value={settings.y?.scale}
            onChange={(event) => { handleChange('y.scale', event.target.value) }}
          >
            {Object.keys(scalesLimited).map((key) =>
              <MenuItem value={key} key={key}>{key}</MenuItem>
            )}
          </TextField>
        </WidgetEditOption>
      </WidgetEditGroup>
      <WidgetEditGroup title="marker color">
        <WidgetEditOption>
          <InputJMESPath
            label="Search quantity"
            value={settings?.markers?.color?.search_quantity}
            onChange={(value) => handleChange('markers.color.search_quantity', value)}
            onSelect={(value) => handleAcceptQuantity('markers.color.search_quantity', value)}
            onAccept={(value) => handleAcceptQuantity('markers.color.search_quantity', value)}
            error={errors['markers.color.search_quantity']}
            onError={(value) => handleErrorQuantity('markers.color.search_quantity', value)}
            dtypes={dtypesColor}
            dtypesRepeatable={dtypesColor}
            optional
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <InputTextField
            label="Title"
            fullWidth
            value={settings.markers?.color?.title}
            onChange={(event) => handleChange('markers.color.title', event.target.value)}
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <UnitInput
            label='Unit'
            value={settings.markers?.color?.unit}
            onChange={(value) => handleChange('markers.color.unit', value)}
            onSelect={(value) => handleAccept('markers.color.unit', value)}
            onAccept={(value) => handleAccept('markers.color.unit', value)}
            error={errors['markers.color.unit']}
            onError={(value) => handleError('markers.color.unit', value)}
            dimension={dimensions['markers.color.search_quantity'] || null}
            optional
            disableGroup
          />
        </WidgetEditOption>
      </WidgetEditGroup>
      <WidgetEditGroup title="General">
        <WidgetEditOption>
          <InputTextField
            label="Title"
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
  x: schemaAxis.required('Search quantity for the x axis is required.'),
  y: schemaAxis.required('Search quantity for the y axis is required.'),
  markers: schemaMarkers,
  size: number().integer().required('Size is required.'),
  autorange: bool()
})
