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
import { Subject } from 'rxjs'
import PropTypes from 'prop-types'
import {
  Box,
  Typography
} from '@material-ui/core'
import { useRecoilValue, RecoilRoot } from 'recoil'
import { convertSI } from '../../utils'
import DOS from './DOS'
import BandStructure from './BandStructure'
import BrillouinZone from './BrillouinZone'
import { unitsState } from '../archive/ArchiveBrowser'
import { makeStyles } from '@material-ui/core/styles'
import { electronicRange } from '../../config'

// Styles
const useStyles = makeStyles((theme) => {
  return {
    row: {
      display: 'flex',
      flexDirection: 'row',
      justifyContent: 'flex-start',
      alignItems: 'flex-start',
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
    },
    noData: {
      top: '1.43rem',
      left: theme.spacing(2),
      right: theme.spacing(2),
      bottom: '3.55rem'
    },
    placeHolder: {
      top: '1.43rem',
      left: theme.spacing(2),
      right: theme.spacing(2),
      bottom: theme.spacing(2)
    }
  }
})

function ElectronicProperties({bs, dos, className, classes}) {
  const units = useRecoilValue(unitsState)
  const range = useMemo(() => convertSI(electronicRange, 'electron_volt', units, false), [units])
  const bsLayout = useMemo(() => ({yaxis: {autorange: false, range: range}}), [range])
  const dosLayout = useMemo(() => ({yaxis: {autorange: false, range: range}}), [range])

  // RxJS subject for efficiently propagating y axis changes between DOS and BS
  const bsYSubject = useMemo(() => new Subject(), [])
  const dosYSubject = useMemo(() => new Subject(), [])

  // Styles
  const styles = useStyles({classes: classes})

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
    <RecoilRoot>
      <Box className={styles.row}>
        <Box className={styles.bs}>
          <Typography variant="subtitle1" align='center'>Band structure</Typography>
          <BandStructure
            data={bs === false ? false : bs?.section_k_band}
            layout={bsLayout}
            aspectRatio={1.2}
            placeHolderStyle={styles.placeHolder}
            noDataStyle={styles.noData}
            unitsState={unitsState}
            onRelayouting={handleBSRelayouting}
            onReset={() => { bsYSubject.next({yaxis: {range: electronicRange}}) }}
            layoutSubject={dosYSubject}
            metaInfoLink={bs?.path}
          ></BandStructure>
        </Box>
        <Box className={styles.dos}>
          <Typography variant="subtitle1" align='center'>Density of states</Typography>
          <DOS
            data={dos === false ? false : dos?.section_dos}
            layout={dosLayout}
            aspectRatio={0.6}
            placeHolderStyle={styles.placeHolder}
            noDataStyle={styles.noData}
            onRelayouting={handleDOSRelayouting}
            onReset={() => { dosYSubject.next({yaxis: {range: electronicRange}}) }}
            unitsState={unitsState}
            layoutSubject={bsYSubject}
            metaInfoLink={dos?.path}
          ></DOS>
        </Box>
        {bs !== false
          ? <Box className={styles.bz}>
            <Typography variant="subtitle1" align='center'>Brillouin zone</Typography>
            <BrillouinZone
              data={bs?.section_k_band}
              aspectRatio={1.2}
            ></BrillouinZone>
          </Box>
          : null
        }
      </Box>
    </RecoilRoot>
  )
}

ElectronicProperties.propTypes = {
  dos: PropTypes.any, // Data for DOS. Set false if no data is available, set to some other falsy value to enable placeholder.
  bs: PropTypes.any, // Data for BS. Set false if no data is available, set to some other falsy value to enable placeholder.
  className: PropTypes.string,
  classes: PropTypes.object
}

export default ElectronicProperties
