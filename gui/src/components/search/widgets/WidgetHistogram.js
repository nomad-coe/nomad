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
import React, { useCallback, useMemo } from 'react'
import PropTypes from 'prop-types'
import { useSearchContext } from '../SearchContext'
import { Widget } from './Widget'
import { ActionCheckbox, ActionSelect } from '../../Actions'
import { Range } from '../input/InputRange'
import { scales } from '../../plotting/common'
import {getDisplayLabel} from '../../../utils'
import { Unit } from '../../units/Unit'
import { useUnitContext } from '../../units/UnitContext'

/**
 * Displays a histogram widget.
 */
export const autorangeDescription = 'Automatically center the view on the data'
export const WidgetHistogram = React.memo((
{
  id,
  title,
  description,
  x,
  nbins,
  scale,
  autorange,
  showinput,
  className
}) => {
  const { filterData, useSetWidget } = useSearchContext()
  const {units} = useUnitContext()
  const setWidget = useSetWidget(id)

  // Create final axis config for the plot
  const xAxis = useMemo(() => {
    const xFilter = filterData[x.quantity]
    const xTitle = x.title || getDisplayLabel(xFilter)
    const xType = xFilter?.dtype
    const xUnit = x.unit
      ? new Unit(x.unit)
      : new Unit(xFilter.unit || 'dimensionless').toSystem(units)

    return {
      ...x,
      title: xTitle,
      unit: xUnit,
      dtype: xType
    }
  }, [filterData, x, units])

  const handleEdit = useCallback(() => {
    setWidget(old => { return {...old, editing: true } })
  }, [setWidget])

  const handleChangeScale = useCallback((value) => {
    setWidget(old => { return {...old, scale: value} })
  }, [setWidget])

  return <Widget
    id={id}
    title={title || 'Histogram'}
    description={description}
    onEdit={handleEdit}
    className={className}
    actions={<>
      <ActionCheckbox
        tooltip={autorangeDescription}
        label="autorange"
        value={autorange}
        onChange={(value) => setWidget(old => ({...old, autorange: value}))}
      />
      <ActionSelect
        value={scale}
        options={Object.keys(scales)}
        tooltip="Statistics scaling"
        onChange={handleChangeScale}
      />
    </>}
  >
    <Range
      xAxis={xAxis}
      visible={true}
      nBins={nbins}
      scale={scale}
      anchored={true}
      autorange={autorange}
      showinput={showinput}
      aggId={id}
    />
  </Widget>
})

WidgetHistogram.propTypes = {
  id: PropTypes.string.isRequired,
  title: PropTypes.string,
  description: PropTypes.string,
  x: PropTypes.object,
  nbins: PropTypes.number,
  scale: PropTypes.string,
  autorange: PropTypes.bool,
  showinput: PropTypes.bool,
  className: PropTypes.string
}
