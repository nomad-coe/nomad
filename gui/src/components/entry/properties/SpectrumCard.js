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
import React from 'react'
import PropTypes from 'prop-types'
import PropertyCard from './PropertyCard'
import { useUnits } from '../../../units'
import EELS from '../../visualization/EELS'
import { resolveRef } from '../../archive/metainfo'

export default function SpectrumCard({entryMetadata, archive}) {
  const units = useUnits()

  const hasData = entryMetadata?.results?.properties?.available_properties?.includes('spectra.eels')

  if (!hasData) {
    return null
  }

  const dataRef = archive?.results?.properties?.spectra?.eels
  const data = dataRef && resolveRef(dataRef, archive)

  return <PropertyCard title="Electron Energy Loss Spectrum">
    <EELS
      data={data}
      layout={{yaxis: {autorange: true}}}
      aspectRatio={2}
      units={units}
    />
  </PropertyCard>
}

SpectrumCard.propTypes = {
  entryMetadata: PropTypes.object.isRequired,
  archive: PropTypes.object
}
