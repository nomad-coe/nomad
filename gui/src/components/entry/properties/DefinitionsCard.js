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
import {PropertyCard} from './PropertyCard'
import { useEntryStore } from '../EntryContext'
import { Box, IconButton, Typography, Link } from '@material-ui/core'
import MoreIcon from '@material-ui/icons/MoreVert'
import { ArchiveButton } from '../../nav/Routes'

const DefinitionsCard = React.memo(function DefinitionsCard({index, archive}) {
  const {entryId, uploadId} = useEntryStore()

  if (!index?.quantities?.includes('definitions')) {
    return null
  }

  const actions = (
    <ArchiveButton
      component={IconButton}
      entryId={entryId} uploadId={uploadId}
      path={'definitions'}
    >
      <MoreIcon/>
    </ArchiveButton>
  )
  const name = archive ? archive.definitions?.name || '' : '...'
  const sections = archive ? archive.definitions?.section_definitions?.length : '...'
  return <PropertyCard title="Schema definitions" action={actions}>
    <Box margin={2}>
      <Box marginBottom={1}>
        <Typography>
          This entry contains a schema package {name} with {sections} section definitions.
        </Typography>
      </Box>
      {archive?.definitions?.section_definitions?.map(sectionDef => {
        const path = 'definitions/' + sectionDef.name
        return <Typography key={sectionDef.name}>
          <ArchiveButton component={Link} path={path} entryId={entryId} uploadId={uploadId}>
            {sectionDef.more?.label || sectionDef.name}
          </ArchiveButton>
        </Typography>
      })}

    </Box>
  </PropertyCard>
})
DefinitionsCard.propTypes = {
  index: PropTypes.object,
  archive: PropTypes.object
}

export default DefinitionsCard
