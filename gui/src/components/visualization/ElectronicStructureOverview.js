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
import { Subject } from 'rxjs'
import PropTypes from 'prop-types'
import {
  Box,
  Typography
} from '@material-ui/core'
import DOS from './DOS'
import NoData from './NoData'
import BandStructure from './BandStructure'
import BrillouinZone from './BrillouinZone'
import { RecoilRoot } from 'recoil'
import { unitsState } from '../archive/ArchiveBrowser'
import { makeStyles } from '@material-ui/core/styles'

function ElectronicStructureOverview({data, range, className, classes, raiseError}) {
  const [dosLayout, setDosLayout] = useState({
    autorange: false,
    yaxis: {range: range}
  })
  const [bsLayout, setBsLayout] = useState({
    autorange: false,
    yaxis: {range: range}
  })

  // RxJS subject for efficiently propagating y axis changes between DOS and BS
  const bsYSubject = useMemo(() => new Subject(), [])
  const dosYSubject = useMemo(() => new Subject(), [])

  // Styles
  const useStyles = makeStyles((theme) => {
    return {
      row: {
        display: 'flex',
        flexDirection: 'row',
        justifyContent: 'flex-start',
        alignItems: 'center',
        width: '100%',
        height: '100%',
        flexWrap: 'wrap'
      },
      bz: {
        flex: '0 0 66.6%'
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
    if (data.dos) {
      let update = {yaxis: {range: [event['yaxis.range[0]'], event['yaxis.range[1]']]}}
      bsYSubject.next(update)
    }
  }, [data, bsYSubject])
  const handleDOSRelayouting = useCallback((event) => {
    if (data.bs) {
      let update = {yaxis: {range: [event['yaxis.range[0]'], event['yaxis.range[1]']]}}
      dosYSubject.next(update)
    }
  }, [data, dosYSubject])

  return (
    <RecoilRoot>
      <Box className={style.row}>
        <Box className={style.bs}>
          <Typography variant="subtitle1" align='center'>Band structure</Typography>
          {data.bs
            ? <BandStructure
              data={data?.bs?.section_k_band}
              layout={bsLayout}
              aspectRatio={1.2}
              unitsState={unitsState}
              onRelayouting={handleBSRelayouting}
              onReset={() => { setDosLayout({yaxis: {range: range}}) }}
              layoutSubject={dosYSubject}
            ></BandStructure>
            : <NoData aspectRatio={1.2}/>
          }
        </Box>
        <Box className={style.dos}>
          <Typography variant="subtitle1" align='center'>Density of states</Typography>
          {data.dos
            ? <DOS
              data={data.dos.section_dos}
              layout={dosLayout}
              aspectRatio={0.6}
              onRelayouting={handleDOSRelayouting}
              onReset={() => { setBsLayout({yaxis: {range: range}}) }}
              unitsState={unitsState}
              layoutSubject={bsYSubject}
            ></DOS>
            : <NoData aspectRatio={0.6}/>
          }
        </Box>
        {data.bs
          ? <Box className={style.bz}>
            <Typography variant="subtitle1" align='center'>Brillouin zone</Typography>
            <BrillouinZone
              data={data.bs.section_k_band}
              aspectRatio={1.2}
            ></BrillouinZone>
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
