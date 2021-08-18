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
import _ from 'lodash'
import PropertyCard from './PropertyCard'
import { toMateriaStructure } from '../../../utils'
import Quantity from '../../Quantity'
import { Box, Grid, makeStyles, Typography } from '@material-ui/core'
import { normalizeDisplayValue, encyclopediaEnabled, appBase } from '../../../config'
import Actions from '../../Actions'
import Structure from '../../visualization/Structure'

const useStyles = makeStyles(theme => ({
  properties: {
    height: '100%',
    display: 'flex',
    justifyContent: 'space-between',
    flexDirection: 'column'
  }
}))

export default function MaterialCard({entryMetadata, archive}) {
  const classes = useStyles()
  const archiveUrl = `/entry/id/${entryMetadata.upload_id}/${entryMetadata.entry_id}/archive`
  const structuresUrl = `${archiveUrl}/results/properties/structures`
  const materialId = entryMetadata.results?.material?.material_id

  let structures = null
  if (archive) {
    structures = []
    const structureKinds = ['original', 'conventional', 'primitive']
    const archiveStructures = archive.results.properties?.structures
    if (archiveStructures) {
      structureKinds.forEach(structureKind => {
        const key = `structure_${structureKind}`
        const archiveStructure = archiveStructures[key]
        const structureUrl = `${structuresUrl}/{key}`
        if (archiveStructure) {
          const materiaStructure = toMateriaStructure(archiveStructure, structureKind, structureUrl)
          structures.push(materiaStructure)
        }
      })
    }
    if (structures.length === 0) {
      structures = false
    }
  }

  return <PropertyCard title="Material">
    <Grid container spacing={1}>
      <Grid item xs={5}>
        <Box className={classes.properties}>
          <Quantity column>
            <Quantity quantity="results.material.chemical_formula_hill" label='formula' noWrap data={entryMetadata}/>
            <Quantity quantity="results.material.structural_type" label='structural type' noWrap data={entryMetadata}/>
            <Quantity quantity="results.material.material_name" label='material name' noWrap data={entryMetadata}/>
            {entryMetadata.results?.material?.symmetry &&
              <Quantity row>
                <Quantity
                  quantity="results.material.symmetry.crystal_system"
                  label='crystal system'
                  noWrap
                  data={entryMetadata}
                />
                <Quantity
                  description="Space group symbol and number"
                  label="space group"
                  noWrap
                  data={entryMetadata}
                >
                  <Typography noWrap>
                    {normalizeDisplayValue(_.get(entryMetadata, 'results.material.symmetry.space_group_symbol'))} ({normalizeDisplayValue(_.get(entryMetadata, 'results.material.symmetry.space_group_number'))})
                  </Typography>
                </Quantity>
              </Quantity>
            }
          </Quantity>
          {encyclopediaEnabled && materialId &&
            <Actions
              justifyContent='flex-start'
              color='primary'
              variant='text'
              size='medium'
              actions={[{
                tooltip: 'View this material in the Encyclopedia',
                content: 'Encyclopedia',
                href: `${appBase}/encyclopedia/#/material/${materialId}`
              }]}
            >
            </Actions>
          }
        </Box>
      </Grid>
      <Grid item xs={7} style={{marginTop: '-2rem'}}>
        <Structure
          data={structures}
          materialType={entryMetadata.results?.material?.structural_type}
          aspectRatio={1.5}
          data-testid="viewer-material"
        />
      </Grid>
    </Grid>
  </PropertyCard>
}

MaterialCard.propTypes = {
  entryMetadata: PropTypes.object.isRequired,
  archive: PropTypes.object
}
