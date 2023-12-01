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
import React, {useState, useEffect, useMemo, useCallback} from 'react'
import PropTypes from 'prop-types'
import { useTheme } from '@material-ui/core/styles'
import { Box, Checkbox, Menu, MenuItem, FormControlLabel } from '@material-ui/core'
import { MoreVert } from '@material-ui/icons'
import Plot from '../plotting/Plot'
import { mergeObjects } from '../../utils'
import { getLineStyles } from '../plotting/common'
import { Quantity } from '../units/Quantity'
import { Unit } from '../units/Unit'
import { useUnitContext } from '../units/UnitContext'
import { withErrorHandler } from '../ErrorHandler'
import { Action } from '../Actions'

const energyUnit = new Unit('joule')

/**
 * Graph for Spectra data.
 */
export const spectraError = 'Could not load Spectra.'
const Spectra = React.memo(({
  data,
  layout,
  className,
  'data-testid': testID,
  ...other
}) => {
  const [finalData, setFinalData] = useState(!data ? data : undefined)
  const theme = useTheme()
  const {units} = useUnitContext()
  const [anchorEl, setAnchorEl] = React.useState(null)
  const [spectraNormalize, setSpectraNormalize] = useState(true)

  // Merge custom layout with default layout
  const tmpLayout = useMemo(() => {
    if (data === undefined) {
      return
    }
    const defaultLayout = {
      showlegend: true,
      legend: {
        x: 1,
        y: 1,
        xanchor: 'right',
        yanchor: 'top'
      },
      yaxis: {
        title: {
          text: 'Intensities (a.u.)'
        }
      },
      xaxis: {
        title: {
          text: `Excitation energies (${energyUnit.toSystem(units).label()})`
        }
      }
    }
    return mergeObjects(layout, defaultLayout)
  }, [layout, units, data])

  // Get the maximum intensity to normalize in the checkbox
  const maximumIntensity = useMemo(() => {
    if (!data) return undefined
    const maxes = data.map(trace => Math.max(...trace.intensities))
    return Math.max(...maxes)
  }, [data])

  // The plotted data is loaded only after the first render as a side effect to
  // avoid freezing the UI
  useEffect(() => {
    if (!data) {
      setFinalData(data)
      return
    }

    // Normalizing spectra
    const plotData = []
    const lineStyles = getLineStyles(data.length, theme)
    if (maximumIntensity !== undefined) {
      for (let i = 0; i < data.length; ++i) {
        const trace = data[i]
        const energies = new Quantity(trace.energies, energyUnit).toSystem(units).value()
        const intensities = spectraNormalize
          ? trace.intensities.map((intensity) => intensity / maximumIntensity)
          : trace.intensities
        plotData.push(
          {
            x: energies,
            y: intensities,
            name: `${data[i].label} ${i}`,
            type: 'scatter',
            mode: 'lines',
            line: lineStyles[i]
          }
        )
      }
    }
    setFinalData(plotData)
  }, [data, theme, units, spectraNormalize, maximumIntensity])

  const open = Boolean(anchorEl)
  const openMenu = useCallback((event) => {
    setAnchorEl(event.currentTarget)
  }, [])
  const closeMenu = useCallback(() => {
    setAnchorEl(null)
  }, [])

  return <Box display="flex" flexDirection="column" height="100%" width="100%">
    <Box flex="1 1 auto">
      <Plot
      data={finalData}
      layout={tmpLayout}
      floatTitle={"Spectra"}
      className={className}
      data-testid={`${testID}`}
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
              setSpectraNormalize(!spectraNormalize)
            }}
            color="primary"
            checked={spectraNormalize}
          />
          }
          label='Normalize intensities'
        />
      </MenuItem>
      </Menu>
    </Box>
  </Box>
})

Spectra.propTypes = {
  data: PropTypes.arrayOf(PropTypes.shape({
    type: PropTypes.string,
    label: PropTypes.string,
    intensities: PropTypes.arrayOf(PropTypes.number),
    energies: PropTypes.arrayOf(PropTypes.number)
  })),
  layout: PropTypes.object,
  className: PropTypes.string,
  'data-testid': PropTypes.string
}
Spectra.defaultProps = {
  'data-testid': 'spectra'
}

export default withErrorHandler(spectraError)(Spectra)
