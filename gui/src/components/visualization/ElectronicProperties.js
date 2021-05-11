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
import { Box } from '@material-ui/core'
import { convertSI } from '../../utils'
import DOS from './DOS'
import BandStructure from './BandStructure'
import BrillouinZone from './BrillouinZone'
import SectionTable from './SectionTable'
import { makeStyles } from '@material-ui/core/styles'
import { electronicRange } from '../../config'
import PropertyContainer from './PropertyContainer'

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
      marginTop: theme.spacing(1.5),
      flex: '0 0 65%'
    },
    gaps: {
      marginTop: theme.spacing(1.5),
      flex: '0 0 35%'
    },
    bs: {
      flex: '0 0 65%'
    },
    dos: {
      flex: '0 0 35%'
    },
    nodata: {
      top: theme.spacing(1),
      left: theme.spacing(2),
      right: theme.spacing(2),
      bottom: '3.55rem'
    },
    placeholder: {
      top: theme.spacing(1),
      left: theme.spacing(2),
      right: theme.spacing(2),
      bottom: theme.spacing(2)
    }
  }
})

// Band gap quantities to show. Saved as const object to prevent re-renders
const bandGapQuantities = {
  index: {label: 'Ch.'},
  band_gap: {label: 'Value'},
  band_gap_type: {label: 'Type', placeholder: 'no gap'}
}

const ElectronicProperties = React.memo(({bs, dos, className, classes, units}) => {
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
    <Box className={styles.row}>
      <PropertyContainer title="Band structure" className={styles.bs}>
        <BandStructure
          data={bs}
          layout={bsLayout}
          aspectRatio={0.6 * 65 / 35}
          placeHolderStyle={styles.placeholder}
          noDataStyle={styles.nodata}
          units={units}
          onRelayouting={handleBSRelayouting}
          onReset={() => { bsYSubject.next({yaxis: {range: electronicRange}}) }}
          layoutSubject={dosYSubject}
          data-testid="bs-electronic"
        ></BandStructure>
      </PropertyContainer>
      <PropertyContainer title="Density of states" className={styles.dos}>
        <DOS
          data={dos}
          layout={dosLayout}
          aspectRatio={0.6}
          placeHolderStyle={styles.placeholder}
          noDataStyle={styles.nodata}
          onRelayouting={handleDOSRelayouting}
          onReset={() => { dosYSubject.next({yaxis: {range: electronicRange}}) }}
          units={units}
          layoutSubject={bsYSubject}
          data-testid="dos-electronic"
        ></DOS>
      </PropertyContainer>
      {bs !== false
        ? <><PropertyContainer title="Brillouin zone" className={styles.bz}>
          <BrillouinZone
            data={bs}
            aspectRatio={0.6 * 65 / 35}
            data-testid="bz-electronic"
          ></BrillouinZone>
        </PropertyContainer>
        <PropertyContainer title="Band gaps" className={styles.gaps}>
          <SectionTable
            horizontal
            section="results.properties.electronic.band_structure_electronic.channel_info"
            quantities={bandGapQuantities}
            data={bs?.channel_info}
            aspectRatio={0.6}
            units={units}
          />
        </PropertyContainer>
        </>
        : null
      }
    </Box>
  )
})

ElectronicProperties.propTypes = {
  dos: PropTypes.any, // Set to false if not available, set to other falsy value to show placeholder.
  bs: PropTypes.any, // Set to false if not available, set to other falsy value to show placeholder.
  className: PropTypes.string,
  classes: PropTypes.object,
  units: PropTypes.object // Contains the unit configuration
}

export default ElectronicProperties
