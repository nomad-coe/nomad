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
import { Quantity, useUnits } from '../../units'
import DOS from './DOS'
import BandStructure from './BandStructure'
import BrillouinZone from './BrillouinZone'
import BandGap from './BandGap'
import { makeStyles } from '@material-ui/core/styles'
import { electronicRange } from '../../config'
import { PropertyGrid, PropertyItem } from '../entry/properties/PropertyCard'

// Styles
const useStyles = makeStyles((theme) => {
  return {
    nodata: {
      top: theme.spacing(0.7),
      bottom: '3.55rem'
    },
    placeholder: {
      top: theme.spacing(0.7),
      bottom: theme.spacing(2)
    }
  }
})

const ElectronicProperties = React.memo(({
  bs,
  dos,
  brillouin_zone,
  band_gap,
  classes
}) => {
  const units = useUnits()
  const range = useMemo(() => new Quantity(electronicRange, 'electron_volt').toSystem(units).value(), [units])
  const bsLayout = useMemo(() => ({yaxis: {autorange: false, range: range}}), [range])
  const dosLayout = useMemo(() => ({yaxis: {autorange: false, range: range}}), [range])

  // RxJS subject for efficiently propagating y axis changes between DOS and BS
  const bsYSubject = useMemo(() => new Subject(), [])
  const dosYSubject = useMemo(() => new Subject(), [])

  // Styles
  const styles = useStyles({classes: classes})

  // Synchronize panning between BS/DOS plots
  const handleBSRelayouting = useCallback((event) => {
    const update = {yaxis: {range: [event['yaxis.range[0]'], event['yaxis.range[1]']]}}
    bsYSubject.next(update)
  }, [bsYSubject])
  const handleDOSRelayouting = useCallback((event) => {
    const update = {yaxis: {range: [event['yaxis.range[0]'], event['yaxis.range[1]']]}}
    dosYSubject.next(update)
  }, [dosYSubject])

  return <PropertyGrid>
    <PropertyItem title="Band structure" xs={8}>
      <BandStructure
        data={bs}
        layout={bsLayout}
        units={units}
        onRelayouting={handleBSRelayouting}
        onReset={() => { bsYSubject.next({yaxis: {range: electronicRange}}) }}
        layoutSubject={dosYSubject}
        data-testid="bs-electronic"
      />
    </PropertyItem>
    <PropertyItem title="Density of states" xs={4}>
      <DOS
        data={dos}
        layout={dosLayout}
        placeHolderStyle={styles.placeholder}
        noDataStyle={styles.nodata}
        onRelayouting={handleDOSRelayouting}
        onReset={() => { dosYSubject.next({yaxis: {range: electronicRange}}) }}
        units={units}
        layoutSubject={bsYSubject}
        data-testid="dos-electronic"
      />
    </PropertyItem>
    <PropertyItem title="Brillouin zone" xs={8}>
      <BrillouinZone
        data={brillouin_zone}
        data-testid="bz-electronic"
      />
    </PropertyItem>
    <PropertyItem title="Band gaps" xs={4}>
      <BandGap data={band_gap}/>
    </PropertyItem>
  </PropertyGrid>
})

ElectronicProperties.propTypes = {
  dos: PropTypes.any,
  bs: PropTypes.any,
  brillouin_zone: PropTypes.any,
  band_gap: PropTypes.any,
  classes: PropTypes.object
}

export default ElectronicProperties
