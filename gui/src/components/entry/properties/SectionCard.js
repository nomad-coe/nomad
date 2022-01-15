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
import React, { useEffect, useMemo } from 'react'
import PropTypes from 'prop-types'
import {PropertyCard} from './PropertyCard'
import SectionEditor from '../../archive/SectionEditor'
import { useEntryContext } from '../EntryContext'
import { Box } from '@material-ui/core'

const SectionCard = React.memo(({sectionDef, getSection, ...props}) => {
  const {archive, requireArchive} = useEntryContext()

  useEffect(() => {
    requireArchive()
  }, [requireArchive])

  const section = useMemo(() => {
    if (archive) {
      return getSection(archive)
    }
  }, [archive, getSection])

  return <PropertyCard title={sectionDef.name}>
    <Box margin={2}>
      <SectionEditor sectionDef={sectionDef} section={section} {...props} />
    </Box>
  </PropertyCard>
})

SectionCard.propTypes = {
  sectionDef: PropTypes.object.isRequired,
  getSection: PropTypes.func.isRequired
}

export default SectionCard
