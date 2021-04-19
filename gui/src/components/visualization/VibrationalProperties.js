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
import { Subject } from 'rxjs'
import {
  Box,
  Typography,
  useTheme
} from '@material-ui/core'
import DOS from './DOS'
import BandStructure from './BandStructure'
import { useRecoilValue, RecoilRoot } from 'recoil'
import { convertSI } from '../../utils'
import { unitsState } from '../archive/ArchiveBrowser'
import { makeStyles } from '@material-ui/core/styles'
import { ErrorHandler } from '../ErrorHandler'
import Plot from './Plot'

const useStyles = makeStyles((theme) => {
  return {
    row: {
      display: 'flex',
      flexDirection: 'row',
      justifyContent: 'flex-start',
      flexWrap: 'wrap',
      alignItems: 'center',
      width: '100%',
      height: '100%'
    },
    free_energy: {
      flex: '0 1 50%'
    },
    heat_capacity: {
      flex: '0 1 50%'
    },
    bs: {
      flex: '0 0 66.6%'
    },
    dos: {
      flex: '0 0 33.3%'
    }
  }
})

export default function VibrationalProperties({bs, dos, freeEnergy, heatCapacity, className, classes, raiseError}) {
  // Find minimum and maximum from DOS/BS. Use this range for both plots.
  const units = useRecoilValue(unitsState)
  const range = useMemo(() => {
    let min
    let max
    const energies = dos?.section_dos?.dos_energies
    if (energies) {
      min = Math.min(...energies)
      max = Math.max(...energies)
    }
    return convertSI([min, max], 'joule', units, false)
  }, [dos, units])

  // RxJS subject for efficiently propagating y axis changes between DOS and BS
  const bsYSubject = useMemo(() => new Subject(), [])
  const dosYSubject = useMemo(() => new Subject(), [])
  const bsLayout = useMemo(() => ({yaxis: {autorange: false, range: range, zeroline: true}}), [range])
  const dosLayout = useMemo(() => ({yaxis: {autorange: false, range: range, zeroline: true}}), [range])

  // Styles
  const styles = useStyles(classes)

  // Synchronize panning between BS/DOS plots
  const handleDOSRelayouting = useCallback((event) => {
    let update = {yaxis: {range: [event['yaxis.range[0]'], event['yaxis.range[1]']]}}
    dosYSubject.next(update)
  }, [dosYSubject])
  const handleBSRelayouting = useCallback((event) => {
    let update = {yaxis: {range: [event['yaxis.range[0]'], event['yaxis.range[1]']]}}
    bsYSubject.next(update)
  }, [bsYSubject])

  const theme = useTheme()
  const heatCapacityData = useMemo(() => {
    if (!heatCapacity) {
      return
    }
    return [{
      x: heatCapacity.temperature,
      y: heatCapacity.thermodynamical_property_heat_capacity_C_v,
      type: 'scatter',
      mode: 'lines',
      line: {
        color: theme.palette.primary.main,
        width: 2
      }
    }]
  }, [heatCapacity, theme])

  const heatCapacityLayout = useMemo(() => {
    return {
      xaxis: {
        title: {
          text: 'Temperature (K)'
        },
        zeroline: false
      },
      yaxis: {
        title: {
          text: 'Heat capacity (J/K)'
        },
        zeroline: false
      }
    }
  }, [])
  const freeEnergyData = useMemo(() => {
    if (!freeEnergy) {
      return
    }
    return [{
      x: freeEnergy.temperature,
      y: freeEnergy.vibrational_free_energy_at_constant_volume,
      type: 'scatter',
      mode: 'lines',
      line: {
        color: theme.palette.primary.main,
        width: 2
      }
    }]
  }, [freeEnergy, theme])

  const freeEnergyLayout = useMemo(() => {
    return {
      xaxis: {
        title: {
          text: 'Temperature (K)'
        },
        zeroline: false
      },
      yaxis: {
        title: {
          text: 'Helmholtz free energy (J)'
        },
        zeroline: false
      }
    }
  }, [])

  return (
    <RecoilRoot>
      <Box className={styles.row}>
        <Box className={styles.bs}>
          <Typography variant="subtitle1" align='center'>Phonon dispersion</Typography>
          <BandStructure
            data={bs === false ? false : bs?.section_k_band}
            layout={bsLayout}
            aspectRatio={1.2}
            unitsState={unitsState}
            onRelayouting={handleBSRelayouting}
            onReset={() => { bsYSubject.next({yaxis: {range: range}}) }}
            layoutSubject={dosYSubject}
            metaInfoLink={bs?.path}
            type='vibrational'
          ></BandStructure>
        </Box>
        <Box className={styles.dos}>
          <Typography variant="subtitle1" align='center'>Phonon density of states</Typography>
          <DOS
            data={dos === false ? false : dos?.section_dos}
            layout={dosLayout}
            aspectRatio={0.6}
            unitsState={unitsState}
            onRelayouting={handleDOSRelayouting}
            onReset={() => { dosYSubject.next({yaxis: {range: range}}) }}
            layoutSubject={bsYSubject}
            metaInfoLink={dos?.path}
            type='vibrational'
          ></DOS>
        </Box>
        <Box className={styles.heat_capacity}>
          <Typography variant="subtitle1" align='center'>Heat capacity</Typography>
          <ErrorHandler message='Could not load heat capacity.'>
            <Plot
              data={heatCapacity && heatCapacityData}
              layout={heatCapacityLayout}
              aspectRatio={1}
              floatTitle="Heat capacity"
              metaInfoLink={heatCapacity?.path}
            >
            </Plot>
          </ErrorHandler>
        </Box>
        <Box className={styles.free_energy}>
          <Typography variant="subtitle1" align='center'>Helmholtz free energy</Typography>
          <ErrorHandler message='Could not load Helmholtz free energy.'>
            <Plot
              data={freeEnergy && freeEnergyData}
              layout={freeEnergyLayout}
              aspectRatio={1}
              floatTitle="Helmholtz free energy"
              metaInfoLink={freeEnergy?.path}
            >
            </Plot>
          </ErrorHandler>
        </Box>
      </Box>
    </RecoilRoot>
  )
}

VibrationalProperties.propTypes = {
  dos: PropTypes.any, // Set to false if not available, set to other falsy value to show placeholder.
  bs: PropTypes.any, // Set to false if not available, set to other falsy value to show placeholder.
  freeEnergy: PropTypes.any, // Set to false if not available, set to other falsy value to show placeholder.
  heatCapacity: PropTypes.any, // Set to false if not available, set to other falsy value to show placeholder.
  range: PropTypes.array,
  className: PropTypes.string,
  classes: PropTypes.object,
  raiseError: PropTypes.func
}
