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
import React, { useState, useCallback } from 'react'
import PropTypes from 'prop-types'
import {
  Box,
  Typography
} from '@material-ui/core'
import DOS from './DOS'
import BandStructure from './BandStructure'
import BrillouinZone from './BrillouinZone'
import Placeholder from '../visualization/Placeholder'
import { RecoilRoot } from 'recoil'
import { unitsState } from '../archive/ArchiveBrowser'
import { makeStyles } from '@material-ui/core/styles'
import { ErrorHandler } from '../ErrorHandler'

function ElectronicStructureOverview({data, range, className, classes, raiseError}) {
  const [dosLayout, setDosLayout] = useState({
    yaxis: {range: range}
  })
  const [bsLayout, setBsLayout] = useState({
    yaxis: {range: range}
  })

  // Styles
  const useStyles = makeStyles((theme) => {
    return {
      row: {
        display: 'flex',
        flexDirection: 'row',
        justifyContent: 'center',
        alignItems: 'center',
        width: '100%',
        height: '100%'
      },
      bz: {
        flex: data.dos ? '0 0 25%' : '0 0 33.3%'
      },
      bs: {
        flex: data.dos ? '0 0 50%' : '0 0 66.6%'
      },
      dos: {
        flex: data.bs ? '0 0 25%' : '0 0 90%'
      }
    }
  })
  const style = useStyles(classes)

  // Synchronize panning between BS/DOS plots
  const handleBSRelayouting = useCallback((event) => {
    let update = {
      yaxis: {
        range: [event['yaxis.range[0]'], event['yaxis.range[1]']]
      }
    }
    setDosLayout(update)
  }, [])
  const handleDOSRelayouting = useCallback((event) => {
    let update = {
      yaxis: {
        range: [event['yaxis.range[0]'], event['yaxis.range[1]']]
      }
    }
    setBsLayout(update)
  }, [])

  return (
    <RecoilRoot>
      <Box className={style.row}>
        {data.bs
          ? <Box className={style.bz}>
            <Typography variant="subtitle1" align='center'>Brillouin zone</Typography>
            {data?.bs?.section_k_band
              ? <ErrorHandler message='Could not load Brilllouin zone.'>
                <BrillouinZone
                  data={data.bs.section_k_band}
                  aspectRatio={0.5}
                ></BrillouinZone>
              </ErrorHandler>
              : <Placeholder className={null} aspectRatio={1.1} variant="rect"></Placeholder>
            }
          </Box>
          : null
        }
        {data.bs
          ? <Box className={style.bs}>
            <Typography variant="subtitle1" align='center'>Band structure</Typography>
            {data?.bs?.section_k_band
              ? <ErrorHandler message='Could not load band structure.'>
                <BandStructure
                  data={data.bs.section_k_band}
                  layout={bsLayout}
                  aspectRatio={1.0}
                  unitsState={unitsState}
                  onRelayouting={handleBSRelayouting}
                ></BandStructure>
              </ErrorHandler>
              : <Placeholder className={null} aspectRatio={1.1} variant="rect"></Placeholder>
            }
          </Box>
          : null
        }
        {data.dos
          ? <Box className={style.dos}>
            <Typography variant="subtitle1" align='center'>Density of states</Typography>
            {data?.dos?.section_dos
              ? <ErrorHandler message='Could not load density of states.'>
                <DOS
                  data={data.dos.section_dos}
                  layout={dosLayout}
                  aspectRatio={0.5}
                  onRelayouting={handleDOSRelayouting}
                  unitsState={unitsState}
                ></DOS>
              </ErrorHandler>
              : <Placeholder className={null} aspectRatio={1.1} variant="rect"></Placeholder>
            }
          </Box>
          : null
        }
      </Box>
    </RecoilRoot>
  )
}

ElectronicStructureOverview.propTypes = {
  data: PropTypes.object,
  range: PropTypes.array,
  className: PropTypes.string,
  classes: PropTypes.object,
  raiseError: PropTypes.func
}
ElectronicStructureOverview.defaultProps = {
  range: [-10, 20]
}

export default ElectronicStructureOverview
