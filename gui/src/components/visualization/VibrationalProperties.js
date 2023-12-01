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
import DOS from './DOS'
import BandStructure from './BandStructure'
import HeatCapacity from './HeatCapacity'
import { PropertyGrid, PropertyItem } from '../entry/properties/PropertyCard'
import HelmholtzFreeEnergy from './HelmholtzFreeEnergy'
import { Quantity } from '../units/Quantity'

const VibrationalProperties = React.memo(({
  bs,
  dos,
  freeEnergy,
  heatCapacity,
  units
}) => {
  // Find minimum and maximum from DOS/BS. Use this range for both plots.
  const range = useMemo(() => {
    let range = [undefined, undefined]
    if (dos) {
      const energies = []
      dos.forEach(d => {
        if (d.energies) energies.push(...d.energies)
      })
      const min = Math.min(...energies)
      const max = Math.max(...energies)
      range = new Quantity([min, max], 'joule').toSystem(units).value()
    }
    return range
  }, [dos, units])

  // RxJS subject for efficiently propagating y axis changes between DOS and BS
  const bsYSubject = useMemo(() => new Subject(), [])
  const dosYSubject = useMemo(() => new Subject(), [])
  const bsLayout = useMemo(() => ({yaxis: {autorange: false, range: range}}), [range])
  const dosLayout = useMemo(() => ({yaxis: {autorange: false, range: range}}), [range])

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
    <PropertyItem title="Phonon dispersion" xs={8}>
      <BandStructure
        data={bs}
        layout={bsLayout}
        units={units}
        onRelayouting={handleBSRelayouting}
        onReset={() => { bsYSubject.next({yaxis: {range: range}}) }}
        layoutSubject={dosYSubject}
        metaInfoLink={bs ? bs[0].m_path : null}
        type="vibrational"
        data-testid="bs-phonon"
      />
    </PropertyItem>
    <PropertyItem title="Phonon density of states" xs={4}>
      <DOS
        data={dos}
        layout={dosLayout}
        onRelayouting={handleDOSRelayouting}
        onReset={() => { dosYSubject.next({yaxis: {range: range}}) }}
        units={units}
        layoutSubject={bsYSubject}
        metaInfoLink={dos ? dos[0].m_path : null}
        type="vibrational"
        data-testid="dos-phonon"
      />
    </PropertyItem>
    <PropertyItem title="Heat capacity" xs={6}>
      <HeatCapacity
        data={heatCapacity}
        units={units}
        data-testid="heat-capacity"
      />
    </PropertyItem>
    <PropertyItem title="Helmholtz free energy" xs={6}>
      <HelmholtzFreeEnergy
        data={freeEnergy}
        units={units}
        data-testid="energy-free"
      />
    </PropertyItem>
  </PropertyGrid>
})

VibrationalProperties.propTypes = {
  dos: PropTypes.any, // Set to false if not available, set to other falsy value to show placeholder.
  bs: PropTypes.any, // Set to false if not available, set to other falsy value to show placeholder.
  freeEnergy: PropTypes.any, // Set to false if not available, set to other falsy value to show placeholder.
  heatCapacity: PropTypes.any, // Set to false if not available, set to other falsy value to show placeholder.
  range: PropTypes.array,
  units: PropTypes.object // Contains the unit configuration
}

export default VibrationalProperties
