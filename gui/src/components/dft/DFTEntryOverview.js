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
import React, { useContext, useState, useEffect, useMemo } from 'react'
import PropTypes from 'prop-types'
import { IconButton, Tooltip, Box, Card, CardContent, Grid, CardHeader, Typography, Link, makeStyles, useTheme } from '@material-ui/core'
import { ToggleButton, ToggleButtonGroup } from '@material-ui/lab'
import ArrowForwardIcon from '@material-ui/icons/ArrowForward'
import { apiContext } from '../api'
import ElectronicStructureOverview from '../visualization/ElectronicStructureOverview'
import ApiDialogButton from '../ApiDialogButton'
import Structure from '../visualization/Structure'
import Plot from '../visualization/Plot'
import Placeholder from '../visualization/Placeholder'
import Quantity from '../Quantity'
import { Link as RouterLink } from 'react-router-dom'
import { DOI } from '../search/DatasetList'
import { domains } from '../domains'
import { errorContext } from '../errors'
import { authorList, convertSI } from '../../utils'
import CollapsibleCard from '../CollapsibleCard'
import { resolveRef } from '../archive/metainfo'
import _, { initial } from 'lodash'

import {appBase, encyclopediaEnabled, normalizeDisplayValue} from '../../config'

const useStyles = makeStyles(theme => ({
  root: {
    maxWidth: '1240px',
    minWidth: '1024px'
  },
  error: {
    marginTop: theme.spacing(2)
  },
  toggle: {
    marginBottom: theme.spacing(1)
  },
  structure: {
    marginTop: theme.spacing(1.5),
    width: '100%',
    height: '16.5rem'
  }
}))

/**
 * Shows an informative overview about the selected entry.
 */
export default function DFTEntryOverview({repo, uploadId, calcId}) {
  const classes = useStyles()
  const {api} = useContext(apiContext)
  const {raiseError} = useContext(errorContext)
  const [electronicStructure, setElectronicStructure] = useState(null)
  const [geoOpt, setGeoOpt] = useState(null)
  const [shownSystem, setShownSystem] = useState('original')
  const [structures, setStructures] = useState(null)
  const materialType = repo?.encyclopedia?.material?.material_type

  // When loaded for the first time, start downloading the archive. Once
  // finished, determine the final layout based on it's contents.TODO: When we
  // have more information stored in the ES index, it can be used to select
  // which parts of the Archive should be downloaded to reduce bandwidth.
  useEffect(() => {
    api.archive(uploadId, calcId).then(data => {
      let structs = new Map()

      // Figure out what properties are present by looping over the SCCS. This
      // information will eventually be directly available in the ES index.
      let dos = null
      let bs = null
      const section_run = data.section_run
      const sccs = section_run[0].section_single_configuration_calculation
      for (let i = sccs.length - 1; i > -1; --i) {
        const scc = sccs[i]
        if (!dos && scc.section_dos) {
          dos = {
            'section_system': scc.single_configuration_calculation_to_system_ref,
            'section_method': scc.single_configuration_calculation_to_system_ref,
            'section_dos': scc.section_dos[scc.section_dos.length - 1]
          }
        }
        if (!bs && scc.section_k_band) {
          bs = {
            'section_system': scc.single_configuration_calculation_to_system_ref,
            'section_method': scc.single_configuration_calculation_to_system_ref,
            'section_k_band': scc.section_k_band[scc.section_k_band.length - 1]
          }
        }
      }
      if (dos || bs) {
        setElectronicStructure({
          'dos': dos, 'bs': bs
        })
      }

      // See if there are workflow results
      const section_wf = data.section_workflow
      if (section_wf) {
        const wfType = section_wf.workflow_type

        // Gather energies from geometry optimization
        if (wfType === 'geometry_optimization') {
          const calculations = section_wf.calculations_ref
          let energies = []
          let initialEnergy = null
          calculations.forEach((ref, i) => {
            const calc = resolveRef(ref, data)
            const e = calc.energy_total
            if (i === 0) {
              initialEnergy = e
            }
            energies.push(e - initialEnergy)
          })
          energies = convertSI(energies, 'joule', {energy: 'electron_volt'}, false)
          setGeoOpt({energies: energies})
        }
      }

      // Get the representative system by looping over systems
      let reprSys = null
      const systems = section_run[0].section_system
      for (let i = systems.length - 1; i > -1; --i) {
        const sys = systems[i]
        if (!reprSys && sys.is_representative) {
          const reprSys = {
            'species': sys.atom_species,
            'cell': sys.lattice_vectors ? convertSI(sys.lattice_vectors, 'meter', {length: 'angstrom'}, false) : undefined,
            'positions': convertSI(sys.atom_positions, 'meter', {length: 'angstrom'}, false),
            'pbc': sys.configuration_periodic_dimensions
          }
          structs.set('original', reprSys)
          break
        }
      }

      // Get the conventional (=normalized) system, if present
      let idealSys = data?.section_metadata?.encyclopedia?.material?.idealized_structure
      if (idealSys) {
        const ideal = {
          'species': idealSys.atom_labels,
          'cell': idealSys.lattice_vectors ? convertSI(idealSys.lattice_vectors, 'meter', {length: 'angstrom'}, false) : undefined,
          'positions': idealSys.atom_positions,
          'fractional': true,
          'pbc': idealSys.periodicity
        }
        structs.set('conventional', ideal)
      }

      setStructures(structs)
    }).catch(error => {
      if (error.name === 'DoesNotExist') {
      } else {
        raiseError(error)
      }
    })
  }, [api, raiseError, uploadId, calcId, setElectronicStructure, setStructures])

  const calcData = repo || {uploadId: uploadId, calcId: calcId}
  const loadingRepo = !repo
  const quantityProps = {data: calcData, loading: loadingRepo}
  const authors = loadingRepo ? null : calcData.authors
  const domain = calcData.domain && domains[calcData.domain]
  const structureToggles = useMemo(() => {
    if (structures) {
      const toggles = []
      structures.forEach((value, key) => {
        toggles.push(<ToggleButton key={key} value={key} aria-label={key}>{key}</ToggleButton>)
      })
      return toggles
    }
    return null
  }, [structures])

  let eSize
  if (electronicStructure) {
    if (electronicStructure.bs && electronicStructure.dos) {
      eSize = 12
    } else if (electronicStructure.bs) {
      eSize = 9
    } else if (electronicStructure.dos) {
      eSize = 4
    }
  }
  const theme = useTheme()

  return (
    <Grid container spacing={2}>
      <Grid item xs={4}>
        <CollapsibleCard
          height={'33rem'}
          title='Material'
          action={encyclopediaEnabled && repo?.encyclopedia?.material?.material_id
            ? <Tooltip title="Show the material of this entry in the NOMAD Encyclopedia.">
              <IconButton href={`${appBase}/encyclopedia/#/material/${repo.encyclopedia.material.material_id}`}><ArrowForwardIcon/></IconButton>
            </Tooltip>
            : null
          }
          content={
            <>
              <Quantity column>
                <Quantity row>
                  <Quantity quantity="formula" label='formula' noWrap data={repo}/>
                  <Quantity quantity="dft.system" label='material type' noWrap data={repo}/>
                  <Quantity quantity="encyclopedia.material.material_name" label='material name' noWrap data={repo}/>
                </Quantity>
                {materialType === 'bulk'
                  ? <Quantity row>
                    <Quantity quantity="dft.crystal_system" label='crystal system' noWrap data={repo}/>
                    <Quantity quantity="dft.spacegroup_symbol" label="spacegroup" noWrap data={repo}>
                      <Typography noWrap>
                        {normalizeDisplayValue(_.get(repo, 'dft.spacegroup_symbol'))} ({normalizeDisplayValue(_.get(repo, 'dft.spacegroup'))})
                      </Typography>
                    </Quantity>
                  </Quantity>
                  : null
                }
              </Quantity>
            </>
          }
          fixedContent={structures
            ? <Box className={classes.structure}>
              <ToggleButtonGroup className={classes.toggle} size="small" exclusive value={shownSystem} onChange={(event, value) => { setShownSystem(value) }} aria-label="text formatting">
                {structureToggles}
              </ToggleButtonGroup>
              <Structure system={structures.get(shownSystem)} aspectRatio={4 / 3} options={{view: {fitMargin: 0.75}}}></Structure>
            </Box>
            : <Placeholder className={classes.structure} variant="rect"></Placeholder>
          }
        ></CollapsibleCard>
      </Grid>

      <Grid item xs={4}>
        <CollapsibleCard
          height={'33rem'}
          title='Method'
          content={
            <>
              <Quantity row>
                <Quantity quantity="dft.code_name" label='code name' noWrap data={repo}/>
                <Quantity quantity="dft.code_version" label='code version' noWrap data={repo}/>
              </Quantity>
              <Quantity row>
                <Quantity quantity="dft.electronic_structure_method" label='electronic structure method' noWrap data={repo}/>
              </Quantity>
              <Quantity row>
                <Quantity quantity="dft.restricted" label='restricted' noWrap data={repo}/>
                <Quantity quantity="dft.closed_shell" label='closed shell' noWrap data={repo}/>
              </Quantity>
              <Quantity row>
                <Quantity quantity="dft.xc_functional" label='xc functional' noWrap data={repo}/>
              </Quantity>
              <Quantity row>
                <Quantity quantity="dft.basis_set" label='basis set type' noWrap data={repo}/>
                <Quantity quantity="dft.cutoff" label='plane wave cutoff' noWrap data={repo}/>
              </Quantity>
              <Quantity row>
                <Quantity quantity="dft.pseudopotential" label='pseudopotential' noWrap data={repo}/>
              </Quantity>
              <Quantity row>
                <Quantity quantity="dft.vdw_method" label='vdw method' noWrap data={repo}/>
              </Quantity>
              <Quantity row>
                <Quantity quantity="dft.relativistic" label='relativistic' noWrap data={repo}/>
              </Quantity>
            </>
          }
        ></CollapsibleCard>
      </Grid>

      <Grid item xs={4}>
        <CollapsibleCard
          height={'33rem'}
          title={'Entry'}
          action={<ApiDialogButton title="Repository JSON" data={calcData} />}
          content={
            <>
              <Quantity column>
                <Quantity quantity='comment' placeholder='no comment' {...quantityProps} />
                <Quantity quantity='references' placeholder='no references' {...quantityProps}>
                  {calcData.references &&
                  <div style={{display: 'inline-grid'}}>
                    {calcData.references.map(ref => <Typography key={ref} noWrap>
                      <Link href={ref}>{ref}</Link>
                    </Typography>)}
                  </div>}
                </Quantity>
                <Quantity quantity='authors' {...quantityProps}>
                  <Typography>
                    {authorList(authors || [])}
                  </Typography>
                </Quantity>
                <Quantity quantity='datasets' placeholder='no datasets' {...quantityProps}>
                  {calcData.datasets &&
                  <div>
                    {calcData.datasets.map(ds => (
                      <Typography key={ds.dataset_id}>
                        <Link component={RouterLink} to={`/dataset/id/${ds.dataset_id}`}>{ds.name}</Link>
                        {ds.doi ? <span>&nbsp; (<DOI doi={ds.doi}/>)</span> : ''}
                      </Typography>))}
                  </div>}
                </Quantity>
              </Quantity>
              <Quantity column style={{maxWidth: 350}}>
                <Quantity quantity="mainfile" noWrap ellipsisFront data={repo} withClipboard />
                <Quantity quantity="calc_id" label={`${domain ? domain.entryLabel : 'entry'} id`} noWrap withClipboard data={repo} />
                <Quantity quantity="encyclopedia.material.material_id" label='material id' noWrap data={repo} withClipboard />
                <Quantity quantity="upload_id" label='upload id' data={repo} noWrap withClipboard />
                <Quantity quantity="upload_time" label='upload time' noWrap data={repo}>
                  <Typography noWrap>
                    {new Date(repo.upload_time).toLocaleString()}
                  </Typography>
                </Quantity>
                <Quantity quantity="raw_id" label='raw id' noWrap hideIfUnavailable data={repo} withClipboard />
                <Quantity quantity="external_id" label='external id' hideIfUnavailable noWrap data={repo} withClipboard />
                <Quantity quantity="last_processing" label='last processing' placeholder="not processed" noWrap data={repo}>
                  <Typography noWrap>
                    {new Date(repo.last_processing).toLocaleString()}
                  </Typography>
                </Quantity>
                <Quantity quantity="last_processing" label='processing version' noWrap placeholder="not processed" data={repo}>
                  <Typography noWrap>
                    {repo.nomad_version}/{repo.nomad_commit}
                  </Typography>
                </Quantity>
              </Quantity>
            </>
          }
        ></CollapsibleCard>
      </Grid>

      {electronicStructure
        ? <Grid item xs={eSize}>
          <Card>
            <CardHeader
              title="Electronic properties"
            />
            <CardContent classes={{root: classes.cardContent}}>
              <Box style={{margin: '0 auto 0 auto', width: '100%', height: '36rem'}}>
                <ElectronicStructureOverview
                  data={electronicStructure}>
                </ElectronicStructureOverview>
              </Box>
            </CardContent>
          </Card>
        </Grid>
        : null
      }
      {geoOpt
        ? <Grid item xs={12}>
          <Card>
            <CardHeader
              title="Geometry optimization"
            />
            <CardContent classes={{root: classes.cardContent}}>
              <Plot
                data={[{
                  x: [...Array(geoOpt.energies.length).keys()],
                  y: geoOpt.energies,
                  type: 'scatter',
                  mode: 'lines',
                  line: {
                    color: theme.palette.primary.main,
                    width: 2
                  }
                }]}
                layout={{
                  xaxis: {
                    title: 'Step number',
                    autorange: true,
                    zeroline: false
                  },
                  yaxis: {
                    title: 'Energy (eV)',
                    autorange: true,
                    zeroline: false
                  }
                }}
                // resetLayout={resetLayout}
                aspectRatio={2}
                floatTitle="Geometry optimizaiton"
              >
              </Plot>
            </CardContent>
          </Card>
        </Grid>
        : null
      }
    </Grid>
  )
}

DFTEntryOverview.propTypes = {
  repo: PropTypes.object.isRequired,
  uploadId: PropTypes.string.isRequired,
  calcId: PropTypes.string.isRequired
}
