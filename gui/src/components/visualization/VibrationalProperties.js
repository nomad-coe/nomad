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
import { Box, Typography } from '@material-ui/core'
import DOS from './DOS'
import BandStructure from './BandStructure'
import HeatCapacity from './HeatCapacity'
import HelmholtzFreeEnergy from './HelmholtzFreeEnergy'
import { convertSI } from '../../utils'
import { makeStyles } from '@material-ui/core/styles'

export default function VibrationalProperties({
  bs,
  dos,
  freeEnergy,
  heatCapacity,
  className,
  classes,
  raiseError,
  units
}) {
  // Find minimum and maximum from DOS/BS. Use this range for both plots.
  const range = useMemo(() => {
    let range = [undefined, undefined]
    if (dos?.energies) {
      const min = Math.min(...dos.energies)
      const max = Math.max(...dos.energies)
      range = convertSI([min, max], 'joule', units, false)
    }
    return range
  }, [dos, units])

  // RxJS subject for efficiently propagating y axis changes between DOS and BS
  const bsYSubject = useMemo(() => new Subject(), [])
  const dosYSubject = useMemo(() => new Subject(), [])
  const bsLayout = useMemo(() => ({yaxis: {autorange: false, range: range}}), [range])
  const dosLayout = useMemo(() => ({yaxis: {autorange: false, range: range}}), [range])

  // Styles
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
  const style = useStyles(classes)

  // Synchronize panning between BS/DOS plots
  const handleBSRelayouting = useCallback((event) => {
    let update = {yaxis: {range: [event['yaxis.range[0]'], event['yaxis.range[1]']]}}
    bsYSubject.next(update)
  }, [bsYSubject])
  const handleDOSRelayouting = useCallback((event) => {
    let update = {yaxis: {range: [event['yaxis.range[0]'], event['yaxis.range[1]']]}}
    dosYSubject.next(update)
  }, [dosYSubject])

  return (
    <Box className={style.row}>
      <Box className={style.bs}>
        <Typography variant="subtitle1" align='center'>Phonon dispersion</Typography>
        <BandStructure
          data={bs}
          layout={bsLayout}
          aspectRatio={1.2}
          units={units}
          onRelayouting={handleBSRelayouting}
          onReset={() => { bsYSubject.next({yaxis: {range: range}}) }}
          layoutSubject={dosYSubject}
          metaInfoLink={bs?.m_path}
          type="vibrational"
          data-testid="bs-phonon"
        ></BandStructure>
      </Box>
      <Box className={style.dos}>
        <Typography variant="subtitle1" align="center">Phonon density of states</Typography>
        <DOS
          data={dos}
          layout={dosLayout}
          aspectRatio={0.6}
          onRelayouting={handleDOSRelayouting}
          onReset={() => { dosYSubject.next({yaxis: {range: range}}) }}
          units={units}
          layoutSubject={bsYSubject}
          metaInfoLink={dos?.m_path}
          type="vibrational"
          data-testid="dos-phonon"
        ></DOS>
      </Box>
      <Box className={style.heat_capacity}>
        <Typography variant="subtitle1" align='center'>Heat capacity</Typography>
        <HeatCapacity
          data={heatCapacity}
          aspectRatio={1}
          units={{...units, 'energy': 'joule'}}
          data-testid="heat-capacity"
        />
      </Box>
      <Box className={style.free_energy}>
        <Typography variant="subtitle1" align='center'>Helmholtz free energy</Typography>
        <HelmholtzFreeEnergy
          data={freeEnergy}
          aspectRatio={1}
          units={{...units, 'energy': 'joule'}}
          data-testid="energy-free"
        />
      </Box>
    </Box>
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
  raiseError: PropTypes.func,
  units: PropTypes.object // Contains the unit configuration
}
