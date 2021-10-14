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
import { get } from 'lodash'
import { PropertyCard } from './PropertyCard'
import { toMateriaStructure } from '../../../utils'
import Quantity from '../../Quantity'
import { Box, Grid, makeStyles, Typography } from '@material-ui/core'
import { normalizeDisplayValue, encyclopediaEnabled } from '../../../config'
import { MaterialButton } from '../../nav/Routes'
import Structure from '../../visualization/Structure'

export function Formula({data}) {
  const formula = (data) => {
    const material = data?.results?.material
    if (!material) {
      return null
    }

    return material.chemical_formula_hill || material.chemical_formula_reduced || material.chemical_formula_descriptive
  }

  return <Quantity
    quantity={formula} label='formula' noWrap data={data}
    description="The chemical formula that describes the simulated system or experiment sample."
  />
}

Formula.propTypes = {
  data: PropTypes.object
}

const useStyles = makeStyles(theme => ({
  properties: {
    height: '100%',
    display: 'flex',
    justifyContent: 'space-between',
    flexDirection: 'column',
    alignItems: 'flex-start'
  }
}))

const MaterialCard = React.memo(({entryMetadata, properties, archive}) => {
  const classes = useStyles()
  const archiveUrl = `/entry/id/${entryMetadata.upload_id}/${entryMetadata.entry_id}/archive`
  const structuresUrl = `${archiveUrl}/results/properties/structures`
  const materialId = entryMetadata.results?.material?.material_id

  let structures = null
  if (archive?.results) {
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
            <Formula data={entryMetadata} />
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
                    {normalizeDisplayValue(get(entryMetadata, 'results.material.symmetry.space_group_symbol'))} ({normalizeDisplayValue(get(entryMetadata, 'results.material.symmetry.space_group_number'))})
                  </Typography>
                </Quantity>
              </Quantity>
            }
          </Quantity>
          {encyclopediaEnabled && materialId &&
            <MaterialButton materialId={materialId}>
              Encyclopedia
            </MaterialButton>
          }
        </Box>
      </Grid>
      <Grid item xs={7} style={{marginTop: '-2rem'}}>
        <Structure
          data={structures}
          materialType={entryMetadata.results?.material?.structural_type}
          aspectRatio={1.4}
          data-testid="viewer-material"
        />
      </Grid>
    </Grid>
  </PropertyCard>
})

MaterialCard.propTypes = {
  entryMetadata: PropTypes.object.isRequired,
  properties: PropTypes.object.isRequired,
  archive: PropTypes.object
}

export default MaterialCard
