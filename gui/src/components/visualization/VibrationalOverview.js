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
import React, { useState, useCallback, useMemo } from 'react'
import PropTypes from 'prop-types'
import {
  Box,
  Typography,
  useTheme
} from '@material-ui/core'
import DOS from './DOS'
import BandStructure from './BandStructure'
import Placeholder from '../visualization/Placeholder'
import { RecoilRoot } from 'recoil'
import { unitsState } from '../archive/ArchiveBrowser'
import { makeStyles } from '@material-ui/core/styles'
import { ErrorHandler } from '../ErrorHandler'
import Plot from '../visualization/Plot'

function VibrationalOverview({data, range, className, classes, raiseError}) {
  const [dosLayout, setDosLayout] = useState({
    yaxis: {
      autorange: true
    }
  })
  const [bsLayout, setBsLayout] = useState({
    yaxis: {
      autorange: true
    }
  })

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
    let update = {
      yaxis: {
        autorange: false,
        range: [event['yaxis.range[0]'], event['yaxis.range[1]']]
      }
    }
    setDosLayout(update)
  }, [])
  const handleDOSRelayouting = useCallback((event) => {
    let update = {
      yaxis: {
        autorange: false,
        range: [event['yaxis.range[0]'], event['yaxis.range[1]']]
      }
    }
    setBsLayout(update)
  }, [])

  const theme = useTheme()
  const heatCapacityData = useMemo(() => {
    return [{
      x: data.temperature,
      y: data.heat_capacity,
      type: 'scatter',
      mode: 'lines',
      line: {
        color: theme.palette.primary.main,
        width: 2
      }
    }]
  }, [data, theme])

  const heatCapacityLayout = useMemo(() => {
    return {
      xaxis: {
        title: 'Temperature (K)',
        zeroline: false
      },
      yaxis: {
        title: 'Heat capacity (J/K)',
        zeroline: false
      }
    }
  }, [])
  const freeEnergyData = useMemo(() => {
    return [{
      x: data.temperature,
      y: data.free_energy,
      type: 'scatter',
      mode: 'lines',
      line: {
        color: theme.palette.primary.main,
        width: 2
      }
    }]
  }, [data, theme])

  const freeEnergyLayout = useMemo(() => {
    return {
      xaxis: {
        title: 'Temperature (K)',
        zeroline: false
      },
      yaxis: {
        title: 'Helmholtz free energy (J)',
        zeroline: false
      }
    }
  }, [])

  return (
    <RecoilRoot>
      <Box className={style.row}>
        {data.bs
          ? <Box className={style.bs}>
            <Typography variant="subtitle1" align='center'>Phonon dispersion</Typography>
            {data?.bs?.section_k_band
              ? <BandStructure
                data={data.bs.section_k_band}
                layout={bsLayout}
                aspectRatio={1.2}
                unitsState={unitsState}
                onRelayouting={handleBSRelayouting}
                onReset={() => { setDosLayout({yaxis: {autorange: true}}) }}
              ></BandStructure>
              : <Placeholder className={null} aspectRatio={1.1} variant="rect"></Placeholder>
            }
          </Box>
          : null
        }
        {data.dos
          ? <Box className={style.dos}>
            <Typography variant="subtitle1" align='center'>Phonon density of states</Typography>
            {data?.dos?.section_dos
              ? <DOS
                data={data.dos.section_dos}
                layout={dosLayout}
                aspectRatio={0.6}
                onRelayouting={handleDOSRelayouting}
                onReset={() => { setBsLayout({yaxis: {autorange: true}}) }}
                unitsState={unitsState}
              ></DOS>
              : <Placeholder className={null} aspectRatio={1.1} variant="rect"></Placeholder>
            }
          </Box>
          : null
        }
        {data.heat_capacity && data.temperature
          ? <Box className={style.heat_capacity}>
            <Typography variant="subtitle1" align='center'>Heat capacity</Typography>
            <ErrorHandler message='Could not load heat capacity.'>
              <Plot
                data={heatCapacityData}
                layout={heatCapacityLayout}
                aspectRatio={1}
                floatTitle="Heat capacity"
              >
              </Plot>
            </ErrorHandler>
          </Box>
          : null
        }
        {data.free_energy && data.temperature
          ? <Box className={style.free_energy}>
            <Typography variant="subtitle1" align='center'>Helmholtz free energy</Typography>
            <ErrorHandler message='Could not load Helmholtz free energy.'>
              <Plot
                data={freeEnergyData}
                layout={freeEnergyLayout}
                aspectRatio={1}
                floatTitle="Helmholtz free energy"
              >
              </Plot>
            </ErrorHandler>
          </Box>
          : null
        }
      </Box>
    </RecoilRoot>
  )
}

VibrationalOverview.propTypes = {
  data: PropTypes.object,
  range: PropTypes.array,
  className: PropTypes.string,
  classes: PropTypes.object,
  raiseError: PropTypes.func
}
VibrationalOverview.defaultProps = {
  range: [-10, 20]
}

export default VibrationalOverview
