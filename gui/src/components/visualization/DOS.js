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
import React, {useEffect, useState, useMemo, useCallback} from 'react'
import PropTypes from 'prop-types'
import { useTheme } from '@material-ui/core/styles'
import { MoreVert } from '@material-ui/icons'
import Plot from '../plotting/Plot'
import { add, mergeObjects } from '../../utils'
import { Quantity, Unit } from '../../units'
import { withErrorHandler } from '../ErrorHandler'
import { Action } from '../Actions'
import { msgNormalizationWarning } from '../../config'
import { getLineStyles } from '../plotting/common'
import { Checkbox, Menu, MenuItem, FormControlLabel } from '@material-ui/core'

const energyUnit = new Unit('joule')
const valueUnit = new Unit('1/joule')

const DOS = React.memo(({
  data,
  layout,
  className,
  units,
  type,
  'data-testid': testID,
  ...other
}) => {
  // Merge custom layout with default layout
  const initialLayout = useMemo(() => {
    const defaultLayout = {
      yaxis: {
        title: {
          text: `Energy (${energyUnit.toSystem(units).label()})`
        },
        zeroline: type === 'vibrational'
      },
      xaxis: {
        showexponent: 'first',
        zeroline: false,
        autorange: false
      },
      showlegend: true,
      legend: {
        x: 1,
        xanchor: 'right',
        yanchor: 'bottom',
        y: 0
      }
    }
    return mergeObjects(layout, defaultLayout)
  }, [layout, units, type])

  const [finalData, setFinalData] = useState(!data ? data : undefined)
  const [finalLayout, setFinalLayout] = useState(initialLayout)
  const [normalizedToHOE, setNormalizedToHOE] = useState(true)
  const theme = useTheme()
  const [dosNormalize, setDosNormalize] = useState(false)
  const [anchorEl, setAnchorEl] = React.useState(null)

  // Side effect that runs when the data that is displayed should change. By
  // running all this heavy stuff within useEffect (instead of e.g. useMemo),
  // the first render containing the placeholders etc. can be done as fast as
  // possible.
  useEffect(() => {
    if (!data) {
      setFinalData(data)
      return
    }

    // Create the final data that will be plotted.
    const plotData = []
    let normalized
    const mins = []
    const maxes = []
    const dosRenormalized = data.map(d => {
      const densities = dosNormalize
      ? d.densities.map((density, index) => {
        return density.map(density_intensity => density_intensity * d.normalization_factors[index])
      })
      : d.densities
      return {...d, densities: densities}
    })
    const channels = dosRenormalized.map(d => d.densities.length)
    const lines = getLineStyles(channels.reduce((a, b) => a + b, 0), null, channels.length > 0 ? channels[0] : null).values()
    dosRenormalized.forEach(d => {
      // Determine the energy reference.
      let energyHighestOccupied
      if (type === 'vibrational') {
        energyHighestOccupied = 0
        normalized = true
      } else {
        if (d.energy_highest_occupied === undefined) {
          energyHighestOccupied = 0
          normalized = false
        } else {
          energyHighestOccupied = new Quantity(d.energy_highest_occupied, energyUnit).toSystem(units).value()
          normalized = true
        }
      }

      // Convert units and determine range
      const nChannels = d.densities.length
      let energies = new Quantity(d.energies, energyUnit).toSystem(units).value()
      const values1 = new Quantity(d.densities[0], valueUnit).toSystem(units).value()
      let values2
      mins.push(Math.min(...values1))
      maxes.push(Math.max(...values1))
      if (nChannels === 2) {
        values2 = new Quantity(d.densities[1], valueUnit).toSystem(units).value()
        mins.push(Math.min(...values2))
        maxes.push(Math.max(...values2))
      }
      if (energyHighestOccupied !== 0) {
        energies = add(energies, -energyHighestOccupied)
      }

      const line = lines.next().value
      if (nChannels === 2) {
        plotData.push(
          {
            x: values2,
            y: energies,
            type: 'scatter',
            mode: 'lines',
            showlegend: false,
            line: lines.next().value
          }
        )
      }
      plotData.push(
        {
          x: values1,
          y: energies,
          name: d.name,
          type: 'scatter',
          mode: 'lines',
          showlegend: d.name !== undefined,
          line: line
        }
      )
    })

    const range = [Math.min(...mins), Math.max(...maxes)]
    // Normalization line
    if (type !== 'vibrational' && normalizedToHOE) {
      plotData.push({
        x: range,
        y: [0, 0],
        name: 'Highest occupied',
        showlegend: true,
        type: 'line',
        mode: 'lines',
        line: {
          color: '#000',
          width: 1
        }
      })
    }

    // Compute the final layout
    const computedLayout = mergeObjects(
      {
        xaxis: {
          title: {
            text: dosNormalize ? `${valueUnit.toSystem(units).label()}` : `states ${valueUnit.toSystem(units).label()}`
          },
          range: range
        },
        legend: {
          font: {
            size: 13
          }
        }
      },
      initialLayout
    )
    setFinalData(plotData)
    setFinalLayout(computedLayout)
    setNormalizedToHOE(normalized)
  }, [data, units, initialLayout, normalizedToHOE, theme, type, dosNormalize])

  const open = Boolean(anchorEl)
  const openMenu = useCallback((event) => {
    setAnchorEl(event.currentTarget)
  }, [])
  const closeMenu = useCallback(() => {
    setAnchorEl(null)
  }, [])

  return <>
    <Plot
      data={finalData}
      layout={finalLayout}
      floatTitle="Density of states"
      warning={normalizedToHOE === false ? msgNormalizationWarning : null}
      metaInfoLink={data ? data[0]?.m_path : null}
      data-testid={testID}
      className={className}
      actions={
        <Action tooltip='Options' onClick={openMenu}>
          <MoreVert/>
        </Action>
      }
      {...other}
    >
    </Plot>
    <Menu
      id='settings-menu'
      anchorEl={anchorEl}
      getContentAnchorEl={null}
      anchorOrigin={{ vertical: 'bottom', horizontal: 'right' }}
      transformOrigin={{ vertical: 'top', horizontal: 'right' }}
      keepMounted
        open={open}
        onClose={closeMenu}
    >
      <MenuItem key='normalization'>
        <FormControlLabel control={
          <Checkbox
            onChange={() => {
              setDosNormalize(data[0].normalization_factors ? (dosNormalize) => !dosNormalize : false)
            }}
            color="primary"
            checked={dosNormalize}
          />
          }
          label='Normalize intensities'
        />
      </MenuItem>
    </Menu>
  </>
})

DOS.propTypes = {
  data: PropTypes.oneOfType([
    PropTypes.bool, // Set to False to show NoData component
    PropTypes.arrayOf(PropTypes.shape({
      energies: PropTypes.array.isRequired, // DOS energies array
      densities: PropTypes.array.isRequired, // DOS values array
      energy_highest_occupied: PropTypes.number, // Highest occupied energy.
      m_path: PropTypes.string, // Path of the section containing the data in the Archive
      label: PropTypes.string // Label of the data
    }))
  ]),
  layout: PropTypes.object,
  className: PropTypes.string,
  units: PropTypes.object, // Contains the unit configuration
  type: PropTypes.string, // Type of band structure: electronic or vibrational
  'data-testid': PropTypes.string
}

DOS.defaultProps = {
  type: 'electronic'
}

export default withErrorHandler('Could not load density of states.')(DOS)
