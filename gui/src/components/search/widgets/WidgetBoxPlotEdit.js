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
import React, { useState, useCallback } from 'react'
import PropTypes from 'prop-types'
import { number, bool, string, reach } from 'yup'
import { cloneDeep } from 'lodash'
import {
  TextField,
  MenuItem,
  Checkbox,
  FormControlLabel
} from '@material-ui/core'
import { useSearchContext } from '../SearchContext'
import { InputJMESPath } from '../input/InputMetainfo'
import { InputTextField } from '../input/InputText'
import UnitInput from '../../units/UnitInput'
import { schemaWidget, schemaAxis, schemaMarkers } from './Widget'
import { WidgetEditDialog, WidgetEditGroup, WidgetEditOption } from './WidgetEdit'
import { DType, parseJMESPath, setDeep, isEmptyString } from '../../../utils'

// Allowed dtypes for different axes and marker color. Predefined in order to not break memoization.
const dtypesNumeric = new Set([DType.Float, DType.Int, DType.Timestamp])
const dtypesCategory = new Set([DType.String, DType.Enum, DType.Int])
const dtypesColor = new Set([DType.String, DType.Enum, DType.Float, DType.Int])

const sampleSizeOptions = {
  100: 100,
  1000: 1000,
  10000: 10000
}
const groupBySizeOptions = {
  5: '5',
  10: '10',
  15: '15',
  20: '20',
  25: '25',
  30: '30'
}

/**
 * A dialog that is used to configure a box plot widget.
 */
export const WidgetBoxPlotEdit = React.memo(({widget}) => {
    const { filterData, useSetWidget } = useSearchContext()
    const [settings, setSettings] = useState(cloneDeep(widget))
    const [errors, setErrors] = useState({})
    const [dimensions, setDimensions] = useState({})
    const setWidget = useSetWidget(widget.id)

    const handleChange = useCallback((key, value) => {
      setSettings(old => {
        const newValue = {...old}
        setDeep(newValue, key, value)
        return newValue
      })
    }, [setSettings])

    const handleError = useCallback((key, value) => {
      setErrors(old => ({...old, [key]: value}))
    }, [setErrors])

    const handleErrorQuantity = useCallback((key, value) => {
      handleError(key, value)
      setDimensions((old) => ({...old, [key]: null}))
    }, [handleError])

    const handleAccept = useCallback((key, value) => {
      try {
        reach(schemaWidgetBoxPlot, key).validateSync(value)
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
      const dimension = filterData[quantity]?.dimension || 'dimensionless'
      setDimensions((old) => ({...old, [key]: dimension}))
    }, [handleAccept, filterData])

    const handleClose = useCallback(() => {
      setWidget(old => ({...old, editing: false}))
    }, [setWidget])

    // Final form validation on accept
    const handleEditAccept = useCallback(() => {
      const independentErrors = Object.values(errors).some(x => !!x)
      if (independentErrors) return

      const yEmpty = isEmptyString(settings?.y?.search_quantity)
      if (yEmpty) {
        handleErrorQuantity('y.search_quantity', 'Please specify a value.')
      }

      if (!independentErrors && !yEmpty) {
        setWidget(old => ({...old, ...{...settings, editing: false, visible: true}}))
      }
    }, [settings, setWidget, errors, handleErrorQuantity])

    return <WidgetEditDialog
        id={widget.id}
        open={widget.editing}
        visible={widget.visible}
        title="Edit box plot widget"
        description='Displays a box plot for a numerical search quantity, optionally grouped and subgrouped by a categorical search quantity.'
        onClose={handleClose}
        onAccept={handleEditAccept}
      >
      <WidgetEditGroup title="Data & display">
        <WidgetEditOption>
          <TextField
            select
            fullWidth
            label="Maximum number of entries to load"
            variant="filled"
            value={settings.sample_size || 100}
            onChange={(event) => { handleChange('sample_size', Number(event.target.value)) }}
          >
            {Object.entries(sampleSizeOptions).map(([key, label]) =>
              <MenuItem value={Number(key)} key={key}>{label}</MenuItem>
            )}
          </TextField>
        </WidgetEditOption>
        <WidgetEditOption>
          <FormControlLabel
            control={<Checkbox checked={settings.show_points} onChange={(event, value) => handleChange('show_points', value)}/>}
            label='Show individual data points.'
          />
          <FormControlLabel
            control={<Checkbox checked={settings.autorange} onChange={(event, value) => handleChange('autorange', value)}/>}
            label='Automatically center the view on the data'
          />
        </WidgetEditOption>
      </WidgetEditGroup>
      <WidgetEditGroup title="y axis">
        <WidgetEditOption>
          <InputJMESPath
            label="Search quantity"
            value={settings.y?.search_quantity}
            error={errors['y.search_quantity']}
            onChange={(value) => handleChange('y.search_quantity', value)}
            onAccept={(value) => handleAcceptQuantity('y.search_quantity', value)}
            onSelect={(value) => handleAcceptQuantity('y.search_quantity', value)}
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
      </WidgetEditGroup>
      <WidgetEditGroup title="Grouping">
        <WidgetEditOption>
          <InputJMESPath
            label="Group by (optional)"
            value={settings.group_by}
            error={errors['group_by']}
            onChange={(value) => handleChange('group_by', value)}
            onAccept={(value) => handleAccept('group_by', value)}
            onSelect={(value) => handleAccept('group_by', value)}
            onError={(value) => handleError('group_by', value)}
            dtypes={dtypesCategory}
            dtypesRepeatable={dtypesCategory}
            optional
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <TextField
            select
            fullWidth
            label="Maximum groups"
            variant="filled"
            value={settings.group_by_size || 10}
            onChange={(event) => { handleChange('group_by_size', event.target.value) }}
          >
            {Object.entries(groupBySizeOptions).map(([key, label]) =>
              <MenuItem value={Number(key)} key={key}>{label}</MenuItem>
            )}
          </TextField>
        </WidgetEditOption>
        <WidgetEditOption>
          <InputJMESPath
            label="Subgroup by (optional)"
            value={settings.subgroup_by}
            error={errors['subgroup_by']}
            onChange={(value) => handleChange('subgroup_by', value)}
            onAccept={(value) => handleAccept('subgroup_by', value)}
            onSelect={(value) => handleAccept('subgroup_by', value)}
            onError={(value) => handleError('subgroup_by', value)}
            dtypes={dtypesCategory}
            dtypesRepeatable={dtypesCategory}
            optional
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <TextField
            select
            fullWidth
            label="Maximum subgroups"
            variant="filled"
            value={settings.subgroup_by_size || 10}
            onChange={(event) => { handleChange('subgroup_by_size', event.target.value) }}
          >
            {Object.entries(groupBySizeOptions).map(([key, label]) =>
              <MenuItem value={Number(key)} key={key}>{label}</MenuItem>
            )}
          </TextField>
        </WidgetEditOption>
      </WidgetEditGroup>
      {settings?.show_points && (
        <WidgetEditGroup title="Marker color">
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
      )}
      <WidgetEditGroup title="Other">
        <WidgetEditOption>
          <InputTextField
            label="Title"
            fullWidth
            value={settings?.title}
            onChange={(event) => handleChange('title', event.target.value)}
          />
        </WidgetEditOption>
        <WidgetEditOption>
          <InputTextField
            label="Description"
            fullWidth
            value={settings?.description}
            multiline
            maxRows={10}
            onChange={(event) => handleChange('description', event.target.value)}
          />
        </WidgetEditOption>
      </WidgetEditGroup>
    </WidgetEditDialog>
})

WidgetBoxPlotEdit.propTypes = {
  widget: PropTypes.object
}

export const schemaWidgetBoxPlot = schemaWidget.shape({
  y: schemaAxis.required('Search quantity for the y axis is required.'),
  sample_size: number().integer().nullable(),
  show_points: bool().nullable(),
  markers: schemaMarkers,
  group_by: string().nullable(),
  group_by_size: number().integer().nullable(),
  subgroup_by: string().nullable(),
  subgroup_by_size: number().integer().nullable(),
  autorange: bool()
})
