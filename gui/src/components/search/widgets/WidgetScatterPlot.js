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
import { string, number, bool } from 'yup'
import { isEmpty } from 'lodash'
import {
  Divider,
  TextField,
  MenuItem,
  Tooltip,
  makeStyles,
  Checkbox,
  FormControlLabel
} from '@material-ui/core'
import { ToggleButton, ToggleButtonGroup } from '@material-ui/lab'
import { InputSearchMetainfo } from '../input/InputText'
import { Widget, schemaWidget } from './Widget'
import { WidgetEditDialog, WidgetEditGroup, WidgetEditOption } from './WidgetEdit'
import { useSearchContext } from '../SearchContext'
import Floatable from '../../visualization/Floatable'
import PlotScatter from '../../plotting/PlotScatter'
import { Action, ActionCheckbox } from '../../Actions'
import { CropFree, PanTool, Fullscreen, Replay } from '@material-ui/icons'
import { autorangeDescription } from './WidgetHistogram'
import { styled } from '@material-ui/core/styles'
import { DType } from '../../../utils'
import { Quantity, Unit, useUnits } from '../../../units'

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
  }
}))

export const WidgetScatterPlot = React.memo((
{
  id,
  label,
  description,
  x,
  y,
  color,
  size,
  autorange,
  dragmode,
  className,
  onSelected
}) => {
  const styles = useStyles()
  const units = useUnits()
  const canvas = useRef()
  const [float, setFloat] = useState(false)
  const [loading, setLoading] = useState(true)
  const { useSetWidget, useHits, filterData, useSetFilter } = useSearchContext()
  const setXFilter = useSetFilter(x)
  const setYFilter = useSetFilter(y)
  const discrete = useMemo(() => {
    return new Set([DType.String, DType.Enum]).has(filterData[color]?.dtype)
  }, [filterData, color])
  const unitX = filterData[x]?.unit || 'dimensionless'
  const unitY = filterData[y]?.unit || 'dimensionless'
  const unitColor = filterData[color]?.unit
  const setWidget = useSetWidget(id)
  const pagination = useMemo(() => ({
    page_size: size,
    order: 'asc',
    order_by: 'entry_id'
  }), [size])
  const required = useMemo(() => {
    const include = ['entry_id']
    !isEmpty(x) && include.push(x)
    !isEmpty(y) && include.push(y)
    !isEmpty(color) && include.push(color)
    return {include}
  }, [x, y, color])

  useEffect(() => {
    setLoading(true)
  }, [required, pagination])

  const hitsCallback = useCallback(() => {
    setLoading(false)
  }, [])
  const hits = useHits(id, required, pagination, hitsCallback)

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
    if (range) {
      const unitXConverted = new Unit(unitX).toSystem(units)
      const unitYConverted = new Unit(unitY).toSystem(units)
      setXFilter({
        gt: new Quantity(range.x[0], unitXConverted),
        lt: new Quantity(range.x[1], unitXConverted)
      })
      setYFilter({
        gt: new Quantity(range.y[0], unitYConverted),
        lt: new Quantity(range.y[1], unitYConverted)
      })
      onSelected && onSelected(data)
    }
  }, [onSelected, setXFilter, setYFilter, unitX, unitY, units])

  const handleDeselect = () => {
    onSelected && onSelected(undefined)
  }

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
      <PlotScatter
        data={loading ? undefined : hits}
        x={x}
        y={y}
        color={color}
        unitX={unitX}
        unitY={unitY}
        unitColor={unitColor}
        discrete={discrete}
        autorange={autorange}
        onSelected={handleSelected}
        onDeselect={handleDeselect}
        dragmode={dragmode}
        onNavigateToEntry={handleNavigated}
        data-testid={id}
        ref={canvas}
      />
    </Widget>
  </Floatable>
})

WidgetScatterPlot.propTypes = {
  id: PropTypes.string.isRequired,
  label: PropTypes.string,
  description: PropTypes.string,
  x: PropTypes.string,
  y: PropTypes.string,
  color: PropTypes.string,
  size: PropTypes.number,
  autorange: PropTypes.bool,
  dragmode: PropTypes.string,
  className: PropTypes.string,
  onSelected: PropTypes.func
}

/**
 * A dialog that is used to configure a scatter plot widget.
 */
export const WidgetScatterPlotEdit = React.memo((props) => {
    const { id, editing, visible } = props
    const { useSetWidget } = useSearchContext()
    const [settings, setSettings] = useState(props)
    const [errors, setErrors] = useState({})
    const setWidget = useSetWidget(id)
    const hasError = useMemo(() => {
      return Object.values(errors).some((d) => !!d) || !schemaWidgetScatterPlot.isValidSync(settings)
    }, [errors, settings])

    const handleSubmit = useCallback((settings) => {
      setWidget(old => ({...old, ...settings}))
    }, [setWidget])

    const handleChange = useCallback((key, value) => {
      setSettings(old => ({...old, [key]: value}))
    }, [setSettings])

    const handleError = useCallback((key, value) => {
      setErrors(old => ({...old, [key]: value}))
    }, [setErrors])

    const handleClose = useCallback(() => {
      setWidget(old => ({...old, editing: false}))
    }, [setWidget])

    const handleAccept = useCallback((key, value) => {
      try {
        schemaWidgetScatterPlot.validateSyncAt(key, {[key]: value})
      } catch (e) {
        handleError(key, e.message)
        return
      }
      setErrors(old => ({...old, [key]: undefined}))
      setSettings(old => ({...old, [key]: value}))
    }, [handleError, setSettings])

    const handleEditAccept = useCallback(() => {
      handleSubmit({...settings, editing: false, visible: true})
    }, [handleSubmit, settings])

    return <WidgetEditDialog
        id={id}
        open={editing}
        visible={visible}
        title="Edit scatter plot widget"
        onClose={handleClose}
        onAccept={handleEditAccept}
        error={hasError}
      >
      <WidgetEditGroup title="x axis">
        <WidgetEditOption>
          <InputSearchMetainfo
            label="quantity"
            value={settings.x}
            error={errors.x}
            onChange={(value) => handleChange('x', value)}
            onSelect={(value) => handleAccept('x', value)}
            onError={(value) => handleError('x', value)}
            dtypes={new Set([DType.Int, DType.Float])}
            repeats={false}
          />
        </WidgetEditOption>
      </WidgetEditGroup>
      <WidgetEditGroup title="y axis">
        <WidgetEditOption>
          <InputSearchMetainfo
            label="quantity"
            value={settings.y}
            error={errors.y}
            onChange={(value) => handleChange('y', value)}
            onSelect={(value) => handleAccept('y', value)}
            onError={(value) => handleError('y', value)}
            dtypes={new Set([DType.Int, DType.Float])}
            repeats={false}
          />
        </WidgetEditOption>
      </WidgetEditGroup>
      <WidgetEditGroup title="color">
        <WidgetEditOption>
          <InputSearchMetainfo
            label="quantity"
            value={settings.color}
            error={errors.color}
            onChange={(value) => handleChange('color', value)}
            onSelect={(value) => handleAccept('color', value)}
            onError={(value) => handleError('color', value)}
            dtypes={new Set([DType.String, DType.Enum, DType.Float, DType.Int])}
            dtypesRepeatable={new Set([DType.String, DType.Enum])}
          />
        </WidgetEditOption>
      </WidgetEditGroup>
      <WidgetEditGroup title="general">
        <WidgetEditOption>
          <TextField
            select
            fullWidth
            label="Maximum number of points"
            variant="filled"
            value={settings.size}
            onChange={(event) => { handleChange('size', event.target.value) }}
          >
            <MenuItem value={100}>100</MenuItem>
            <MenuItem value={1000}>1000</MenuItem>
            <MenuItem value={10000}>10000</MenuItem>
          </TextField>
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
  id: PropTypes.string.isRequired,
  editing: PropTypes.bool,
  visible: PropTypes.bool,
  x: PropTypes.string,
  y: PropTypes.string,
  color: PropTypes.string,
  size: PropTypes.number,
  autorange: PropTypes.bool,
  onClose: PropTypes.func
}

export const schemaWidgetScatterPlot = schemaWidget.shape({
  x: string().required('Quantity for the x axis is required.'),
  y: string().required('Quantity for the y axis is required.'),
  color: string(),
  size: number().integer().required('Size is required.'),
  autorange: bool()
})
