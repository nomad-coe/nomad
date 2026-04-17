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
import { number, bool } from 'yup'
import { isEmpty, range } from 'lodash'
import {
  Divider,
  Tooltip,
  makeStyles
} from '@material-ui/core'
import { ToggleButton, ToggleButtonGroup, Alert } from '@material-ui/lab'
import { Widget, schemaWidget, schemaAxis, schemaMarkers } from './Widget'
import { useSearchContext } from '../SearchContext'
import Floatable from '../../visualization/Floatable'
import PlotScatter from '../../plotting/PlotScatter'
import { Action } from '../../Actions'
import { Search, CropFree, PanTool, Fullscreen, Replay } from '@material-ui/icons'
import { styled } from '@material-ui/core/styles'
import {DType, parseJMESPath} from '../../../utils'
import { Quantity } from '../../units/Quantity'
import { Unit } from '../../units/Unit'
import { useUnitContext } from '../../units/UnitContext'
import { getAxisConfig, getData } from '../../plotting/common'

const stripDtypeSuffix = (quantity) => {
  return quantity?.replace(/#(int|float|str|bool|Datetime)$/, '') || quantity
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
  title,
  description,
  x,
  y,
  markers,
  size,
  autorange,
  drag_mode,
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
  const [xParsed, yParsed, colorParsed, discrete, error] = useMemo(() => {
    const xParsed = parseJMESPath(x?.search_quantity)
    const yParsed = parseJMESPath(y?.search_quantity)
    const colorParsed = markers?.color?.search_quantity ? parseJMESPath(markers.color.search_quantity) : {}
    const discrete = colorParsed?.quantity && new Set([DType.String, DType.Enum]).has(filterData[colorParsed.quantity]?.dtype)
    if (xParsed.error || yParsed.error || colorParsed.error) {
      return [{}, {}, {}, undefined, 'Invalid JMESPath query, please check your syntax.']
    }
    return [xParsed, yParsed, colorParsed, discrete, undefined]
  }, [markers?.color?.search_quantity, x?.search_quantity, y?.search_quantity, filterData])

  // Get storage unit for API communication
  const {storageUnitX, storageUnitY, storageUnitColor} = useMemo(() => {
    if (error) return {}
    const storageUnitX = new Unit(filterData[xParsed.quantity]?.unit || 'dimensionless')
    const storageUnitY = new Unit(filterData[yParsed.quantity]?.unit || 'dimensionless')
    const storageUnitColor = new Unit(filterData[colorParsed?.quantity]?.unit || 'dimensionless')
    return {storageUnitX, storageUnitY, storageUnitColor}
  }, [filterData, xParsed.quantity, yParsed.quantity, colorParsed?.quantity, error])

  // Create final axis config for the plot
  const {xAxis, yAxis, colorAxis} = useMemo(() => {
    if (error) return {xAxis: {}, yAxis: {}, colorAxis: {}}
    return {
      xAxis: getAxisConfig(x, filterData, units),
      yAxis: getAxisConfig(y, filterData, units),
      colorAxis: markers?.color ? getAxisConfig(markers.color, filterData, units) : {}
    }
  }, [error, filterData, markers?.color, x, units, y])

  const setXFilter = useSetFilter(stripDtypeSuffix(xParsed.quantity))
  const setYFilter = useSetFilter(stripDtypeSuffix(yParsed.quantity))

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
    const x = []
    const y = []
    const color = colorParsed.path ? [] : undefined
    const id = []
    for (const hit of hits) {
      const {hitData, error, nPoints} = getData(hit, xParsed.path, yParsed.path, colorParsed.path, discrete)
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
    const x = xAxis.dtype === DType.Timestamp
      ? dataRaw.x
      : new Quantity(dataRaw.x, storageUnitX).to(xAxis.unit).value()
    const y = yAxis.dtype === DType.Timestamp
      ? dataRaw.y
      : new Quantity(dataRaw.y, storageUnitY).to(yAxis.unit).value()
    const color = dataRaw.color && (discrete
      ? dataRaw.color
      : new Quantity(dataRaw.color, storageUnitColor).to(colorAxis.unit).value()
    )
    return {x, y, color, id: dataRaw.id}
  }, [dataRaw, xAxis.dtype, xAxis.unit, storageUnitX, yAxis.dtype, yAxis.unit, storageUnitY, discrete, storageUnitColor, colorAxis.unit])

  const handleEdit = useCallback(() => {
    setWidget(old => { return {...old, editing: true } })
  }, [setWidget])

  const handleDragModeChanged = useCallback((event, value) => {
    if (value !== null) {
      setWidget(old => ({...old, drag_mode: value}))
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
      gte: new Quantity(range.x[0], xAxis.unit),
      lte: new Quantity(range.x[1], xAxis.unit)
    })
    setYFilter({
      gte: new Quantity(range.y[0], yAxis.unit),
      lte: new Quantity(range.y[1], yAxis.unit)
    })
    onSelected?.(data)
  }, [setXFilter, xAxis.unit, setYFilter, yAxis.unit, onSelected])

  const handleDeselect = useCallback(() => {
    onSelected?.(undefined)
  }, [onSelected])

  const actions = useMemo(() => {
      return <>
        <StyledToggleButtonGroup
          size="small"
          value={drag_mode}
          exclusive
          onChange={handleDragModeChanged}
        >
          <ToggleButton value="zoom">
            <Tooltip title="Zoom">
              <Search fontSize="small"/>
            </Tooltip>
          </ToggleButton>
          <ToggleButton value="pan">
            <Tooltip title="Pan">
              <PanTool fontSize="small"/>
            </Tooltip>
          </ToggleButton>
          <ToggleButton value="select">
            <Tooltip title="Select region to create a filter">
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
  }, [drag_mode, handleDragModeChanged, handleResetClick, handleFloat, styles])

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
      title={title || "Scatter plot"}
      description={description}
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
          dragMode={drag_mode}
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
  title: PropTypes.string,
  description: PropTypes.string,
  x: PropTypes.object,
  y: PropTypes.object,
  markers: PropTypes.object,
  size: PropTypes.number,
  autorange: PropTypes.bool,
  drag_mode: PropTypes.string,
  className: PropTypes.string,
  onSelected: PropTypes.func
}

export const schemaWidgetScatterPlot = schemaWidget.shape({
  x: schemaAxis.required('Quantity for the x axis is required.'),
  y: schemaAxis.required('Quantity for the y axis is required.'),
  markers: schemaMarkers,
  size: number().integer().required('Size is required.'),
  autorange: bool()
})
