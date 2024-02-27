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
import React, {useCallback, useEffect, useMemo, useRef, useState} from 'react'
import PropTypes from 'prop-types'
import { number, bool, reach } from 'yup'
import jmespath from 'jmespath'
import { isEmpty, isArray, isEqual, cloneDeep, range, isNil, flattenDeep } from 'lodash'
import {
  Divider,
  Tooltip,
  makeStyles,
  Checkbox,
  FormControlLabel
} from '@material-ui/core'
import { ToggleButton, ToggleButtonGroup, Alert } from '@material-ui/lab'
import { InputJMESPath } from '../input/InputMetainfo'
import { Widget, schemaWidget, schemaAxis, schemaMarkers } from './Widget'
import { WidgetEditDialog, WidgetEditGroup, WidgetEditOption, WidgetEditSelect } from './WidgetEdit'
import { useSearchContext } from '../SearchContext'
import Floatable from '../../visualization/Floatable'
import PlotScatter from '../../plotting/PlotScatter'
import { Action, ActionCheckbox } from '../../Actions'
import { CropFree, PanTool, Fullscreen, Replay } from '@material-ui/icons'
import { autorangeDescription } from './WidgetHistogram'
import { styled } from '@material-ui/core/styles'
import { DType, setDeep, parseJMESPath } from '../../../utils'
import { Quantity } from '../../units/Quantity'
import { Unit } from '../../units/Unit'
import { useUnitContext } from '../../units/UnitContext'
import { InputTextField } from '../input/InputText'
import UnitInput from '../../units/UnitInput'

// Predefined in order to not break memoization
const dtypesNumeric = new Set([DType.Int, DType.Float])
const dtypesColor = new Set([DType.String, DType.Enum, DType.Float, DType.Int])
const nPointsOptions = {
  100: 100,
  1000: 1000,
  10000: 10000
}

const StyledToggleButtonGroup = styled(ToggleButtonGroup)(({ theme }) => ({
  '& .MuiToggleButtonGroup-grouped': {
    margin: theme.spacing(0, 0.25),
    border: 0,
    '&.Mui-disabled': {
      border: 0
    },
    '&:not(:first-of-type)': {
      borderRadius: theme.shape.borderRadius
    },
    '&:first-of-type': {
      borderRadius: theme.shape.borderRadius
    }
  }
}))

/**
 * A thin wrapper for displaying a Plotly scatter plot as a widget.
 */
const useStyles = makeStyles((theme) => ({
  widget: {
    width: '100%',
    height: '100%'
  },
  divider: {
    margin: theme.spacing(0.5, 0.5)
  },
  alert: {
    overflow: 'auto'
  }
}))

export const WidgetScatterPlot = React.memo((
{
  id,
  label,
  description,
  x,
  y,
  markers,
  size,
  autorange,
  dragmode,
  className,
  onSelected
}) => {
  const styles = useStyles()
  const {units} = useUnitContext()
  const canvas = useRef()
  const [float, setFloat] = useState(false)
  const [loading, setLoading] = useState(true)
  const { useSetWidget, useHits, filterData, useSetFilter } = useSearchContext()

  // Parse additional JMESPath config
  const [xParsed, yParsed, colorParsed, error] = useMemo(() => {
    const xParsed = parseJMESPath(x?.quantity)
    const yParsed = parseJMESPath(y?.quantity)
    const colorParsed = markers?.color?.quantity ? parseJMESPath(markers.color.quantity) : {}
    if (xParsed.error || yParsed.error || colorParsed.error) {
      return [{}, {}, {}, 'Invalid JMESPath query, please check your syntax.']
    }
    return [xParsed, yParsed, colorParsed, undefined]
  }, [markers?.color?.quantity, x?.quantity, y?.quantity])

  // Parse units
  const {unitXObj, unitYObj, unitColorObj, displayUnitX, displayUnitY, displayUnitColor, discrete} = useMemo(() => {
    if (error) return {}
    const unitXObj = new Unit(filterData[xParsed.quantity].unit || 'dimensionless')
    const unitYObj = new Unit(filterData[yParsed.quantity].unit || 'dimensionless')
    const unitColorObj = new Unit(filterData[colorParsed?.quantity]?.unit || 'dimensionless')
    const displayUnitX = x.unit ? new Unit(x.unit) : unitXObj.toSystem(units)
    const displayUnitY = y.unit ? new Unit(y.unit) : unitYObj.toSystem(units)
    const displayUnitColor = markers?.color?.unit ? new Unit(markers?.color?.unit) : unitColorObj.toSystem(units)
    const discrete = colorParsed?.quantity && new Set([DType.String, DType.Enum]).has(filterData[colorParsed.quantity]?.dtype)
    return {unitXObj, unitYObj, unitColorObj, displayUnitX, displayUnitY, displayUnitColor, discrete}
  }, [filterData, x.unit, xParsed.quantity, y.unit, yParsed.quantity, markers?.color?.unit, colorParsed?.quantity, units, error])

  // Create final axis config for the plot
  const {xAxis, yAxis, colorAxis} = useMemo(() => {
    if (error) return {}
    const xTitle = x.title || filterData[xParsed.quantity]?.label
    const yTitle = y.title || filterData[yParsed.quantity]?.label
    const colorTitle = markers?.color?.title || filterData[colorParsed.quantity]?.label
    const unitLabelX = displayUnitX.label()
    const unitLabelY = displayUnitY.label()
    const unitLabelColor = displayUnitColor.label()
    return {
      xAxis: {...x, ...xParsed, title: xTitle, unit: unitLabelX},
      yAxis: {...y, ...yParsed, title: yTitle, unit: unitLabelY},
      colorAxis: markers?.color ? {...markers.color, ...colorParsed, title: colorTitle, unit: unitLabelColor} : {}
    }
  }, [colorParsed, displayUnitColor, displayUnitX, displayUnitY, filterData, markers?.color, x, xParsed, y, yParsed, error])

  const setXFilter = useSetFilter(xParsed.quantity)
  const setYFilter = useSetFilter(yParsed.quantity)

  const setWidget = useSetWidget(id)
  const pagination = useMemo(() => ({
    page_size: size,
    order: 'asc'
  }), [size])
  const required = useMemo(() => {
    const include = new Set(['entry_id'])
    for (const config of [xParsed, yParsed, colorParsed]) {
      if (!isEmpty(config)) {
        include.add(config.quantity)
        for (const extra of config.extras) {
          include.add(extra)
        }
      }
    }
    return {include: [...include]}
  }, [xParsed, yParsed, colorParsed])

  useEffect(() => {
    setLoading(true)
  }, [required, pagination])

  const hitsCallback = useCallback(() => {
    setLoading(false)
  }, [])

  // Fetch the data using the useHits hook that automatically applies the
  // existing filters. We filter out data that is invalid (e.g. no values or
  // incompatible sizes between x/y/color). TODO: The API should support an
  // "exists" query that could be used to return hits that actually have the
  // requested values.  This way we would get a better match for the query size
  // and we would not need to manually validate and check the results.
  const hits = useHits(id, required, pagination, hitsCallback)
  const dataRaw = useMemo(() => {
    if (!hits || error) return
    function getData(hit) {
      const hitData = {}

      // Get each property using JMESPath. Errors at this stage will simply
      // cause the entry to be ignored.
      for (const [name, path] of [['x', xParsed.path], ['y', yParsed.path], ['color', colorParsed.path]]) {
        if (isEmpty(path)) continue
        let value
        try {
          value = jmespath.search(hit, path)
        } catch (e) {
          return {error: 'Invalid JMESPATH'}
        }
        // Missing x/y/color value will cause an error unless dealing with
        // discretized colors
        if (isNil(value)) {
          if (name === 'color' && discrete) {
            value = 'undefined'
          } else {
            return {error: 'Empty value'}
          }
        }
        hitData[name] = value
      }

      // Get the shapes
      const xShape = getShape(hitData.x)
      const yShape = getShape(hitData.y)

      // Check if x/y leaf shapes match
      if (xShape[xShape.length - 1] !== yShape[yShape.length - 1]) {
        return {error: 'Incompatible size for x/y'}
      }

      // If x/y shapes do not match, extend accordingly
      const biggestShape = [xShape, yShape].reduce((prev, current) => {
        return (prev.length > current.length)
          ? prev
          : current
      })
      hitData.x = extendFront(hitData.x, xShape, biggestShape)
      hitData.y = extendFront(hitData.y, yShape, biggestShape)

      // Modify color dimensions
      let colorShape = colorParsed.path && getShape(hitData.color)
      if (colorShape && !isEqual(colorShape, biggestShape)) {
        // If color has one more dimension than other arrays and it is discrete,
        // we reduce the last dimension to a single string
        if (discrete && colorShape.length === biggestShape.length + 1) {
          hitData.color = reduceInner(hitData.color)
          colorShape = colorShape.slice(0, -1)
        }
        // Scalar color values are extended
        if (colorShape.length === 0 || (colorShape.length === 1 && colorShape[0] === 1)) {
          hitData.color = fill(
            biggestShape,
            colorShape.length === 0
               ? hitData.color
               : hitData.color[0]
          )
        // Colors are extended according to traces
        } else if ((colorShape.length < biggestShape.length) && colorShape[0] === biggestShape[0]) {
          hitData.color = extendBack(hitData.color, colorShape, biggestShape)
        } else {
          return {error: 'Incompatible size for color'}
        }
      }

      // Flatten arrays
      hitData.x = flatten(hitData.x)
      hitData.y = flatten(hitData.y)
      hitData.color = colorParsed.path && flatten(hitData.color)

      // If shapes still don't match, skip entry. TODO: This check is not ideal,
      // since we may be accepting accidentally mathing sizes. A proper shape
      // check that would also allow "ragged arrays" would be better.
      if (hitData.x.length !== hitData.y.length || (colorParsed.path && hitData.x.length !== hitData.color.length)) {
        return {error: 'Incompatible number of elements'}
      }

      return {hitData, nPoints: hitData.x.length}
    }
    const x = []
    const y = []
    const color = colorParsed.path ? [] : undefined
    const id = []
    for (const hit of hits) {
      const {hitData, error, nPoints} = getData(hit)
      if (error || !nPoints) continue
      for (const i of range(nPoints)) {
        x.push(hitData.x[i])
        y.push(hitData.y[i])
        colorParsed.path && color.push(hitData.color?.[i])
        id.push(hit.entry_id)
      }
    }
    return {x, y, color, id}
  }, [discrete, hits, xParsed.path, yParsed.path, colorParsed.path, error])

  // Perform unit conversion, report errors
  const data = useMemo(() => {
    if (!dataRaw) return
    const x = new Quantity(dataRaw.x, unitXObj).to(displayUnitX).value()
    const y = new Quantity(dataRaw.y, unitYObj).to(displayUnitY).value()
    const color = dataRaw.color && (discrete
      ? dataRaw.color
      : new Quantity(dataRaw.color, unitColorObj).to(displayUnitColor).value()
    )
    return {x, y, color, id: dataRaw.id}
  }, [dataRaw, displayUnitColor, displayUnitX, displayUnitY, unitColorObj, unitXObj, unitYObj, discrete])

  const handleEdit = useCallback(() => {
    setWidget(old => { return {...old, editing: true } })
  }, [setWidget])

  const handleDragModeChanged = useCallback((event, value) => {
    if (value !== null) {
      setWidget(old => ({...old, dragmode: value}))
    }
  }, [setWidget])

  const handleResetClick = useCallback(() => {
    canvas.current.reset()
  }, [])

  const handleFloat = useCallback(() => {
    // The current layout needs to be saved, because the DOM change will
    // cause it to be lost.
    canvas.current?.saveLayout()
    setFloat(old => !old)
  }, [])

  const handleSelected = useCallback((data) => {
    const range = data?.range
    if (!range) return
    setXFilter({
      gte: new Quantity(range.x[0], displayUnitX),
      lte: new Quantity(range.x[1], displayUnitX)
    })
    setYFilter({
      gte: new Quantity(range.y[0], displayUnitY),
      lte: new Quantity(range.y[1], displayUnitY)
    })
    onSelected?.(data)
  }, [onSelected, setXFilter, setYFilter, displayUnitX, displayUnitY])

  const handleDeselect = useCallback(() => {
    onSelected?.(undefined)
  }, [onSelected])

  const actions = useMemo(() => {
      return <>
        <ActionCheckbox
          tooltip={autorangeDescription}
          label="autorange"
          value={autorange}
          onChange={(value) => setWidget(old => ({...old, autorange: value}))}
        />
        <Divider flexItem orientation="vertical" className={styles.divider} />
        <StyledToggleButtonGroup
          size="small"
          value={dragmode}
          exclusive
          onChange={handleDragModeChanged}
        >
          <ToggleButton value="pan">
            <Tooltip title="Pan">
              <PanTool fontSize="small"/>
            </Tooltip>
          </ToggleButton>
          <ToggleButton value="select">
            <Tooltip title="Focus on region">
              <CropFree fontSize="small"/>
            </Tooltip>
          </ToggleButton>
        </StyledToggleButtonGroup>
        <Divider flexItem orientation="vertical" className={styles.divider} />
        <Action tooltip='Reset view' onClick={handleResetClick}>
          <Replay fontSize="small"/>
        </Action>
        <Action tooltip='Toggle fullscreen' onClick={handleFloat}>
          <Fullscreen fontSize="small"/>
        </Action>
      </>
  }, [dragmode, handleDragModeChanged, handleResetClick, handleFloat, autorange, setWidget, styles])

  const handleNavigated = useCallback(() => {
    setFloat(false)
  }, [])

  return <Floatable
      className={className}
      float={float}
      onFloat={handleFloat}
    >
    <Widget
      id={id}
      label={label || "Scatter plot"}
      description={description || 'Custom scatter plot'}
      onEdit={handleEdit}
      actions={actions}
      className={styles.widget}
    >
      {error
        ? <Alert severity="error" className={styles.alert}>
            {error}
          </Alert>
        : <PlotScatter
          data={loading ? undefined : data}
          xAxis={xAxis}
          yAxis={yAxis}
          colorAxis={colorAxis}
          discrete={discrete}
          autorange={autorange}
          onSelected={handleSelected}
          onDeselect={handleDeselect}
          dragmode={dragmode}
          onNavigateToEntry={handleNavigated}
          data-testid={id}
          ref={canvas}
        />
      }
    </Widget>
  </Floatable>
})

WidgetScatterPlot.propTypes = {
  id: PropTypes.string.isRequired,
  label: PropTypes.string,
  description: PropTypes.string,
  x: PropTypes.object,
  y: PropTypes.object,
  markers: PropTypes.object,
  size: PropTypes.number,
  autorange: PropTypes.bool,
  dragmode: PropTypes.string,
  className: PropTypes.string,
  onSelected: PropTypes.func
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
      if (!independentErrors) {
        setWidget(old => ({...old, ...{...settings, editing: false, visible: true}}))
      }
    }, [settings, setWidget, errors])

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
            onSelect={(value) => handleAccept('markers.color.quantity', value)}
            onAccept={(value) => handleAccept('markers.color.quantity', value)}
            error={errors['markers.color.quantity']}
            onError={(value) => handleError('markers.color.quantity', value)}
            dtypes={dtypesColor}
            dtypesRepeatable={dtypesColor}
            optional
          />
        </WidgetEditOption>
      </WidgetEditGroup>
      <WidgetEditGroup title="general">
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

/**
 * Used to flatten the input into a single array of values.
 */
function flatten(input) {
  return isArray(input)
    ? flattenDeep(input)
    : [input]
}

/**
 * Gets the shape of an abitrarily nested array.
 */
function getShape(input) {
  if (!isArray(input)) {
    return []
  }

  const shape = []
  let inner = input

  while (isArray(inner)) {
    shape.push(inner.length)
    inner = inner[0]
  }

  return shape
}

/**
 * Reduces the innermost dimension into a single value.
 */
function reduceInner(input) {
  function reduceRec(inp) {
    if (isArray(inp)) {
      if (isArray(inp[0])) {
        for (let i = 0; i < inp.length; ++i) {
          inp[i] = reduceRec(inp[i])
        }
      } else {
        return inp.sort().join(", ")
      }
    }
    return inp
  }
  return reduceRec(input)
}

/**
 * Resizes the given array to a new size by extending the data to fit the front
 * dimensions (similar to array = array[None, :] in NumPy).
 */
function extendFront(array, oldShape, newShape) {
  // If shape is already correct, return the input array
  const diff = newShape.length - oldShape.length
  if (diff === 0) return array

  // Extend the array
  const extendedArray = []
  function extendRec(depth) {
    const dim = newShape[depth]
    const hasData = depth === diff
    if (hasData) {
      return array
    } else {
      for (let j = 0; j < dim; ++j) {
        extendedArray.push(extendRec(depth + 1))
      }
    }
    return extendedArray
  }
  extendRec(0)
  return extendedArray
}

/**
 * Resizes the given array to a new size by extending the data to fit the last
 * dimensions dimensions (similar to array = array[:, None] in NumPy).
 */
function extendBack(array, oldShape, newShape) {
  // If shape is already correct, return the input array
  const diff = newShape.length - oldShape.length
  if (diff === 0) return array

  // Extend the array
  const extendedArray = []
  const nTraces = newShape[0]
  for (let i = 0; i < nTraces; ++i) {
    const traceArray = []
    for (let j = 0; j < newShape[1]; ++j) {
      traceArray.push([array[i]])
    }
    extendedArray.push(traceArray)
  }
  return extendedArray
}

/**
 * Creates a new array with the given shape, filled with the given value.
 */
function fill(shape, fillValue) {
  if (shape.length === 0) {
    return fillValue
  } else {
    const innerShape = shape.slice(1)
    const innerArray = []
    for (let i = 0; i < shape[0]; i++) {
        innerArray.push(fill(innerShape, fillValue))
    }
    return innerArray
  }
}
