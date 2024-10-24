/* eslint-disable quotes */
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
import { resolveNomadUrlNoThrow } from '../../../utils'
import { useEntryStore } from '../EntryContext'
import { useMetainfoDef } from '../../archive/metainfo'
import { H5WebView } from '../../archive/ArchiveBrowser'
import { PropertyCard } from './PropertyCard'

const HDF5DatasetCard = React.memo(function HDF5DatasetCard({archive}) {
  const {url, uploadId} = useEntryStore()
  const m_def = archive?.data?.m_def_id ? `${archive.data.m_def}@${archive.data.m_def_id}` : archive?.data?.m_def
  const dataMetainfoDefUrl = url && resolveNomadUrlNoThrow(m_def, url)
  const dataMetainfoDef = useMetainfoDef(dataMetainfoDefUrl)

  if (!dataMetainfoDef || !dataMetainfoDef.m_annotations?.h5web) {
    return null
  }

  return (
    <PropertyCard title='HDF5 Data'>
      <H5WebView section={archive.data} def={dataMetainfoDef} uploadId={uploadId}></H5WebView>
    </PropertyCard>
  )
})
HDF5DatasetCard.propTypes = {
  archive: PropTypes.object.isRequired
}

export default HDF5DatasetCard
