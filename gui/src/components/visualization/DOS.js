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
import { add, mergeObjects, resolveInternalRef } from '../../utils'
import { Quantity } from '../units/Quantity'
import { Unit } from '../units/Unit'
import { useUnitContext } from '../units/UnitContext'
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
  type,
  'data-testid': testID,
  ...other
}) => {
  const {units} = useUnitContext()

  // Merge custom layout with default layout
  const initialLayout = useMemo(() => {
    const defaultLayout = {
      xaxis: {
        showexponent: 'first',
        autorange: false,
        zeroline: type === 'vibrational'
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
  }, [layout, type])

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
    const version = data.length > 0 ? data[0].version : undefined
    const processedData = version === 'new'
      ? processNewDosData(data, units, initialLayout, normalizedToHOE, theme, type, dosNormalize)
      : processOldDosData(data, units, initialLayout, normalizedToHOE, theme, type, dosNormalize)

    if (processedData) {
      setFinalData(processedData.finalData)
      setFinalLayout(processedData.finalLayout)
      setNormalizedToHOE(processedData.normalized)
    }
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
        <Action tooltip='Options' onClick={openMenu} data-testid="dos-options-menu">
          <MoreVert/>
        </Action>
      }
      {...other}
    >
    </Plot>
    <Menu
      anchorEl={anchorEl}
      getContentAnchorEl={null}
      anchorOrigin={{ vertical: 'bottom', horizontal: 'right' }}
      transformOrigin={{ vertical: 'top', horizontal: 'right' }}
      keepMounted
      open={open}
      onClose={closeMenu}
    >
      <MenuItem key='normalization'>
        <FormControlLabel
          control={
            <Checkbox
              onChange={() => {
                setDosNormalize(data[0].normalization_factors ? (dosNormalize) => !dosNormalize : false)
              }}
              color="primary"
              checked={dosNormalize}
              disabled={!data?.[0]?.normalization_factors || data[0].normalization_factors.every(value => value === undefined)}
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
  type: PropTypes.string, // Type of band structure: electronic or vibrational
  'data-testid': PropTypes.string
}

DOS.defaultProps = {
  type: 'electronic'
}

export default withErrorHandler('Could not load density of states.')(DOS)

/**
 * Given index and archive, resolves the data for plotting the DOS.
 */
export function resolveDosNew(properties, archive, pattern) {
  // Return undefined if property not available
  if (!properties.has('dos_electronic_new')) return false

  // Resolve the final data once archive is here
  let dosReferences = archive?.results?.properties?.electronic?.dos_electronic_new || []
  if (!Array.isArray(dosReferences)) dosReferences = [dosReferences]
  const dos = []
  for (const reference of dosReferences) {
    if (reference.spin_polarized && reference.data.length !== 2) break
    const d = {
      energies: [],
      densities: [],
      normalization_factors: [],
      energy_highest_occupied: [],
      version: 'new'
    }
    d.name = reference.label
    d.spin_polarized = reference.spin_polarized
    for (const data of reference.data) {
      const match = data.energies.match(pattern)
      const path = match ? match[2] : data.energies
      const totalMatch = data.total.match(pattern)
      const totalPath = match ? totalMatch[2] : data.total
      const sourceArchive = match ? (archive.m_ref_archives[match[1]] || archive.m_ref_archives[data.energies.split('#')[0]]) : archive
      if (sourceArchive) {
        d.energies.push(resolveInternalRef(path, sourceArchive))
        const internalRef = resolveInternalRef(totalPath, sourceArchive)
        d.densities.push(internalRef.value)
        d.normalization_factors.push(internalRef.normalization_factor)
      }
      if (data.energy_ref) {
        d.energy_highest_occupied.push(data.energy_ref)
      }
      d.m_path = `${archive?.metadata?.entry_id}/data/results/properties/electronic/dos_electronic_new`
    }
    if (d.energies && d.densities) dos.push(d)
  }
  return dos.length === 0 ? false : dos
}

/**
 * Given index and archive, resolves the data for plotting the old DOS. This will eventually
 * be deprecated and fully replaced with resolveDosNew().
 */
export function resolveDosOld(properties, archive, pattern) {
  // Return undefined if property not available
  if (!properties.has('dos_electronic')) return false

  // Resolve the final data once archive is here
  let dosReferences = archive?.results?.properties?.electronic?.dos_electronic || []
  if (!Array.isArray(dosReferences)) dosReferences = [dosReferences]
  const dos = []
  for (const reference of dosReferences) {
    const d = {version: 'old'}
    const match = reference.energies.match(pattern)
    const path = match ? match[2] : reference.energies
    const totalPath = match ? reference.total.map(ref => ref.match(pattern)[2]) : reference.total
    const sourceArchive = match ? (archive.m_ref_archives[match[1]] || archive.m_ref_archives[reference.energies.split('#')[0]]) : archive
    if (sourceArchive) {
      d.energies = resolveInternalRef(path, sourceArchive)
      const internalRef = resolveInternalRef(totalPath, sourceArchive)
      d.densities = internalRef.map(dos => dos.value)
      d.normalization_factors = internalRef.map(dos => dos.normalization_factor)
    }
    d.name = reference.label
    if (reference.band_gap) {
      d.energy_highest_occupied = Math.max(...reference.band_gap.map(x => x.energy_highest_occupied))
    }
    d.m_path = `${archive?.metadata?.entry_id}/data/results/properties/electronic/dos_electronic`
    if (d.energies && d.densities) dos.push(d)
  }
  return dos.length === 0 ? false : dos
}

/**
 * Processed the old DOS data plotting.
 */
function processOldDosData(data, units, initialLayout, normalizedToHOE, theme, type, dosNormalize) {
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
      yaxis: {
        title: {
          text: `Energy (${energyUnit.toSystem(units).label()})`
        }
      },
      xaxis: {
        title: {
          text: dosNormalize ? `${valueUnit.toSystem(units).label()}` : `states ${valueUnit.toSystem(units).label()}`
        },
        range: range,
        zeroline: true
      },
      legend: {
        font: {
          size: 13
        }
      }
    },
    initialLayout
  )

  return {
    finalData: plotData,
    finalLayout: computedLayout,
    normalized: normalized
  }
}
/*
 * Processed the new DOS data plotting.
 */
function processNewDosData(data, units, initialLayout, normalizedToHOE, theme, type, dosNormalize) {
  // Plotting styles: lineStyles differentiates spin channels, while lineColors changes
  // for different methodologies
  // @TODO use getLineStyles instead
  const lineStyles = ['solid', 'dot']
  const lineColors = ['#005271', '#f57c00']

  const plotData = []
  let normalized
  const isSpinPolarized = []
  const mins = []
  const maxes = []
  for (const [iMethod, dos_method] of data.entries()) {
    const methodName = dos_method.name ? dos_method.name : ''
    const spinChannels = dos_method.spin_polarized ? 2 : 1
    isSpinPolarized.push(spinChannels === 2)
    const spinName = spinChannels === 2 ? ['+', '-'] : ['']
    // Looping over different spin channels
    for (let iSpin = 0; iSpin < spinChannels; iSpin++) {
      // Setting state for renormalizing densities by the normalization factor
      const densities = dos_method.densities[iSpin]
      const normalization_factor = dos_method.normalization_factors[iSpin]
      let dosRenormalized
      if (normalization_factor !== undefined && dosNormalize) {
        dosRenormalized = densities.map(d => d * normalization_factor)
      } else {
        dosRenormalized = densities
      }

      // Setting the energy reference using the energy_highest_occupied
      let energyHighestOccupied
      if (type === 'vibrational') {
        energyHighestOccupied = 0
        normalized = true
      } else {
        if (dos_method.energy_highest_occupied[iSpin] === undefined) {
          energyHighestOccupied = 0
          normalized = false
        } else {
          energyHighestOccupied = new Quantity(dos_method.energy_highest_occupied[iSpin], energyUnit).toSystem(units).value()
          normalized = true
        }
      }

      // Convert units and determine range
      let energies = new Quantity(dos_method.energies[iSpin], energyUnit).toSystem(units).value()
      if (energyHighestOccupied !== 0) {
        energies = add(energies, -energyHighestOccupied)
      }
      const densitiesUnits = new Quantity(dosRenormalized, valueUnit).toSystem(units).value()
      mins.push(Math.min(...densitiesUnits))
      maxes.push(Math.max(...densitiesUnits))
      const line = {
        width: 2,
        dash: lineStyles[iSpin],
        color: lineColors[iMethod]
      }
      plotData.push(
        {
          x: densitiesUnits,
          y: energies,
          name: `${methodName}${spinName[iSpin]}`,
          type: 'scatter',
          mode: 'lines',
          showlegend: (methodName !== '' || spinChannels === 2),
          line: line
        }
      )
    }
  }

  // Zero energy line in case Fermi energy (EF) is stored
  const range = [Math.min(...mins), Math.max(...maxes)]
  let referencedLayout
  if (normalized) {
    referencedLayout = {
      yaxis: {
        title: {
          text: `Energy - E<sub>F</sub> (${energyUnit.toSystem(units).label()})`
        },
        zeroline: true
      }
    }
  } else {
    referencedLayout = {
      yaxis: {
        title: {
          text: `Energy (${energyUnit.toSystem(units).label()})`
        }
      }
    }
  }
  const zeroEnergyRefLayout = mergeObjects(initialLayout, referencedLayout)

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
    zeroEnergyRefLayout
  )

  return {
    finalData: plotData,
    finalLayout: computedLayout,
    normalized: normalized
  }
}
