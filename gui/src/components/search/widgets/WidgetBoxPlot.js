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
import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import PropTypes from 'prop-types'
import { isEmpty, range } from 'lodash'
import jmespath from 'jmespath'
import { Divider, Tooltip, makeStyles } from '@material-ui/core'
import { ToggleButton, ToggleButtonGroup } from '@material-ui/lab'
import { styled } from '@material-ui/core/styles'
import { Search, PanTool, Replay, Fullscreen } from '@material-ui/icons'
import Floatable from '../../visualization/Floatable'
import { Widget } from './Widget'
import { Action } from '../../Actions'
import { useSearchContext } from '../SearchContext'
import PlotBox from '../../plotting/PlotBox'
import { getAxisConfig, getData } from '../../plotting/common'
import { useUnitContext } from '../../units/UnitContext'
import { Quantity } from '../../units/Quantity'
import { Unit } from '../../units/Unit'
import { DType, parseJMESPath } from '../../../utils'

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

const useStyles = makeStyles((theme) => ({
  widget: {
    width: '100%',
    height: '100%'
  },
  divider: {
    margin: theme.spacing(0.5, 0.5)
  }
}))

export const WidgetBoxPlot = React.memo((
{
  id,
  title,
  description,
  y,
  markers,
  sample_size,
  show_points,
  group_by,
  group_by_size,
  subgroup_by,
  subgroup_by_size,
  autorange,
  drag_mode = 'zoom',
  className
}) => {
  const styles = useStyles()
  const canvas = useRef()
  const { units } = useUnitContext()
  const { filterData, useSetWidget, useHits } = useSearchContext()
  const setWidget = useSetWidget(id)
  const [loading, setLoading] = useState(true)
  const [float, setFloat] = useState(false)

  // Parse additional JMESPath config
  const [xParsed, yParsed, colorParsed, subgroupParsed, discrete, error] = useMemo(() => {
    const yParsed = parseJMESPath(y?.search_quantity)
    const xParsed = group_by ? parseJMESPath(group_by) : {}
    const colorParsed = markers?.color?.search_quantity ? parseJMESPath(markers.color.search_quantity) : {}
    const subgroupParsed = subgroup_by ? parseJMESPath(subgroup_by) : {}
    const discrete = colorParsed?.quantity && new Set([DType.String, DType.Enum]).has(filterData[colorParsed.quantity]?.dtype)
    return [xParsed, yParsed, colorParsed, subgroupParsed, discrete, undefined]
  }, [markers?.color?.search_quantity, group_by, subgroup_by, y?.search_quantity, filterData])

  // Get storage unit for API communication
  const {storageUnitX, storageUnitY, storageUnitColor} = useMemo(() => {
    if (error) return {}
    const storageUnitY = new Unit(filterData[yParsed.quantity]?.unit || 'dimensionless')
    const storageUnitColor = new Unit(filterData[colorParsed?.quantity]?.unit || 'dimensionless')
    return {storageUnitX, storageUnitY, storageUnitColor}
  }, [filterData, yParsed.quantity, colorParsed?.quantity, error])

  // Create final axis configs for the plot
  const {xAxis, yAxis, colorAxis} = useMemo(() => {
    if (error) return {xAxis: {}, yAxis: {}, colorAxis: {}}
    return {
      xAxis: getAxisConfig({search_quantity: group_by}, filterData, units),
      yAxis: getAxisConfig(y, filterData, units),
      colorAxis: markers?.color ? getAxisConfig(markers.color, filterData, units) : {}
    }
  }, [error, filterData, group_by, markers?.color, units, y])

  // Sampling configuration
  const {pagination, required } = useMemo(() => {
    const include = new Set(['entry_id'])
    for (const config of [xParsed, yParsed, colorParsed, subgroupParsed]) {
      if (config && !isEmpty(config)) {
        include.add(config.quantity)
        for (const extra of config.extras) {
          include.add(extra)
        }
      }
    }

    return {
      pagination: {
        page_size: sample_size || 100,
        order: 'asc'
      },
      required: { include: [...include] }
    }
  }, [yParsed, xParsed, colorParsed, subgroupParsed, sample_size])

  useEffect(() => {
    setLoading(true)
  }, [required, pagination])

  const hitsCallback = useCallback(() => {
    setLoading(false)
  }, [])

  // Fetch hits
  const hits = useHits(
    id,
    required,
    pagination,
    hitsCallback
  )

  const handleEdit = useCallback(() => {
    setWidget(old => ({ ...old, editing: true }))
  }, [setWidget])

  const handleResetClick = useCallback(() => {
    canvas.current?.reset()
  }, [])

  const handleFloat = useCallback(() => {
    canvas.current?.saveLayout()
    setFloat(old => !old)
  }, [])

  const handleDragModeChanged = useCallback((event, value) => {
    if (value !== null) {
      setWidget(old => ({...old, drag_mode: value}))
    }
  }, [setWidget])

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

  // Convert a value from storage unit to display unit
  const convert = useCallback((value) => {
    if (value == null || !yAxis.unit) return value
    return new Quantity(value, storageUnitY).to(yAxis.unit).value()
  }, [storageUnitY, yAxis.unit])

  // Resolve the data for the box plot based on the mode
  const boxData = useMemo(() => {
    // Sample mode – build flat lists, convert units, then group by label
    // so that PlotBox receives a uniform list of groups in both modes.
    if (!hits) return null
    if (!yParsed.path) return null
    const x = []
    let y = []
    let color = colorParsed.path ? [] : undefined
    const subgroup = subgroupParsed.path ? [] : undefined
    const id = []
    for (const hit of hits) {
      const {hitData, error, nPoints} = getData(hit, xParsed?.path, yParsed?.path, colorParsed?.path, discrete)
      if (error || !nPoints) continue
      let subgroupValue
      if (subgroupParsed.path) {
        try {
          subgroupValue = jmespath.search(hit, subgroupParsed.path)
        } catch (e) {
          subgroupValue = undefined
        }
        if (subgroupValue == null) continue
        subgroupValue = String(subgroupValue)
      }
      for (const i of range(nPoints)) {
        x.push(hitData.x[i])
        y.push(hitData.y[i])
        colorParsed.path && color.push(hitData.color?.[i])
        if (subgroup) subgroup.push(subgroupValue)
        id.push(hit.entry_id)
      }
    }

    // Perform unit conversion on flat arrays
    y = (yAxis.dtype === DType.Timestamp)
      ? y
      : new Quantity(y, storageUnitY).to(yAxis.unit).value()
    color = color && (discrete
      ? color
      : new Quantity(color, storageUnitColor).to(colorAxis.unit).value()
    )

    // Group the converted data by label (and optionally subgroup label),
    // sort by count descending, and apply size limits.
    const groups = new Map()
    for (let i = 0; i < y.length; i++) {
      const groupLabel = x[i] ?? undefined
      const subLabel = subgroup ? subgroup[i] : undefined
      const key = subLabel != null
        ? `${groupLabel ?? '__default__'}::${subLabel}`
        : (groupLabel ?? '__default__')
      if (!groups.has(key)) {
        groups.set(key, { label: groupLabel, subgroupLabel: subLabel, y: [], color: color ? [] : undefined, id: [] })
      }
      const group = groups.get(key)
      group.y.push(y[i])
      if (color) group.color.push(color[i])
      group.id.push(id[i])
    }

    if (subgroup) {
      // Collect subgroups under their parent groups
      const outerGroups = new Map()
      for (const item of groups.values()) {
        const key = item.label ?? '__default__'
        if (!outerGroups.has(key)) outerGroups.set(key, [])
        outerGroups.get(key).push(item)
      }
      // Sort outer groups by total point count descending
      const sortedOuter = [...outerGroups.entries()]
        .sort((a, b) => {
          const countA = a[1].reduce((sum, g) => sum + g.y.length, 0)
          const countB = b[1].reduce((sum, g) => sum + g.y.length, 0)
          return countB - countA
        })
      const maxGroups = group_by_size || 10
      const maxSub = subgroup_by_size || 10
      const result = []
      for (const [key, subs] of sortedOuter.slice(0, maxGroups)) {
        subs.sort((a, b) => b.y.length - a.y.length)
        result.push({
          label: key === '__default__' ? undefined : key,
          subgroups: subs.slice(0, maxSub).map(s => ({
            label: s.subgroupLabel,
            y: s.y,
            color: s.color,
            id: s.id
          }))
        })
      }
      return result
    } else {
      const sorted = [...groups.values()].sort((a, b) => b.y.length - a.y.length)
      const maxGroups = group_by_size || 10
      return sorted.slice(0, maxGroups)
    }
  }, [group_by, group_by_size, subgroup_by_size, convert, hits, yParsed?.path, colorParsed?.path, subgroupParsed?.path, yAxis.dtype, yAxis.unit, storageUnitY, discrete, storageUnitColor, colorAxis?.unit, xParsed?.path])

  return <Floatable
      className={className}
      float={float}
      onFloat={handleFloat}
    >
    <Widget
      id={id}
      title={title || 'Box plot'}
      description={description}
      onEdit={handleEdit}
      actions={actions}
      className={styles.widget}
    >
      <div className={styles.widget}>
        <PlotBox
          data={loading ? null : boxData}
          title={group_by ? (title || 'Box plot') : undefined}
          xAxis={group_by ? xAxis : undefined}
          yAxis={yAxis}
          colorAxis={colorAxis}
          discrete={discrete}
          show_points={show_points}
          autorange={autorange}
          dragMode={drag_mode}
          ref={canvas}
        />
      </div>
    </Widget>
  </Floatable>
})

WidgetBoxPlot.propTypes = {
  id: PropTypes.string.isRequired,
  title: PropTypes.string,
  description: PropTypes.string,
  y: PropTypes.object,
  markers: PropTypes.object,
  sample_size: PropTypes.number,
  show_points: PropTypes.bool,
  group_by: PropTypes.string,
  group_by_size: PropTypes.number,
  subgroup_by: PropTypes.string,
  subgroup_by_size: PropTypes.number,
  autorange: PropTypes.bool,
  drag_mode: PropTypes.string,
  className: PropTypes.string
}
