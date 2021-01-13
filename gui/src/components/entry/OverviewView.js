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
import React, { useContext, useState, useEffect } from 'react'
import PropTypes from 'prop-types'
import { Skeleton } from '@material-ui/lab'
import { Tooltip, Box, Card, CardContent, Grid, CardHeader, Typography, Link, makeStyles } from '@material-ui/core'
import { apiContext } from '../api'
import ApiDialogButton from '../ApiDialogButton'
import ElectronicStructureOverview from '../visualization/ElectronicStructureOverview'
import Structure from '../visualization/Structure'
import Quantity from '../Quantity'
import { Link as RouterLink } from 'react-router-dom'
import { DOI } from '../search/DatasetList'
import { domains } from '../domains'
import { EntryPageContent } from './EntryPage'
import { errorContext } from '../errors'
import { authorList, convertSI } from '../../utils'
import _ from 'lodash'
import {appBase, encyclopediaEnabled, normalizeDisplayValue} from '../../config'

const useStyles = makeStyles(theme => ({
  root: {
    maxWidth: '1240px',
    minWidth: '1024px'
  },
  error: {
    marginTop: theme.spacing(2)
  },
  cardContent: {
    paddingTop: 0
  },
  topCard: {
    height: '32rem'
  },
  structureViewer: {
    height: '25rem',
    padding: '0'
  }
}))

/**
 * Shows an informative overview about the selected entry.
 */
export default function OverviewView({uploadId, calcId}) {
  const classes = useStyles()
  const {api} = useContext(apiContext)
  const {raiseError} = useContext(errorContext)
  const [repo, setRepo] = useState(null)
  const [exists, setExists] = useState(true)
  const [electronicStructure, setElectronicStructure] = useState(null)
  const [reprSys, setReprSys] = useState(null)

  // When loaded for the first time, download calc data from the ElasticSearch
  // index.
  useEffect(() => {
    api.repo(uploadId, calcId).then(data => {
      console.log('Loaded repo!')
      setRepo(data)
    }).catch(error => {
      if (error.name === 'DoesNotExist') {
        setExists(false)
      } else {
        raiseError(error)
      }
    })
  }, [api, raiseError, uploadId, calcId, setRepo, setExists])

  // When loaded for the first time, start downloading the archive. Once
  // finished, determine the final layout based on it's contents.TODO: When we
  // have more information stored in the ES index, it can be used to select
  // which parts of the Archive should be downloaded to reduce bandwidth.
  useEffect(() => {
    api.archive(uploadId, calcId).then(data => {
      console.log('Loaded archive!')

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
      setElectronicStructure({
        'dos': dos, 'bs': bs
      })

      // Get the representative system by looping over systems
      let reprSys = null
      const systems = section_run[0].section_system
      for (let i = systems.length - 1; i > -1; --i) {
        const sys = systems[i]
        if (!reprSys && sys.is_representative) {
          reprSys = {
            'species': sys.atom_species,
            'cell': sys.lattice_vectors ? convertSI(sys.lattice_vectors, 'meter', {length: 'angstrom'}, false) : undefined,
            'positions': convertSI(sys.atom_positions, 'meter', {length: 'angstrom'}, false),
            'pbc': sys.configuration_periodic_dimensions
          }
          break
        }
      }
      setReprSys(reprSys)
    }).catch(error => {
      if (error.name === 'DoesNotExist') {
        setExists(false)
      } else {
        raiseError(error)
      }
    })
  }, [api, raiseError, uploadId, calcId, setExists, setElectronicStructure, setReprSys])

  const calcData = repo || {uploadId: uploadId, calcId: calcId}
  const loadingRepo = !repo
  const quantityProps = {data: calcData, loading: loadingRepo}

  const authors = loadingRepo ? null : calcData.authors
  const domain = calcData.domain && domains[calcData.domain]

  const {data} = repo || {uploadId: uploadId, calcId: calcId}
  const material_name = entry => {
    let name
    try {
      name = entry.encyclopedia.material.material_name
    } catch {}
    name = name || 'unnamed'

    if (encyclopediaEnabled && data?.encyclopedia?.material?.material_id && data.published && !data.with_embargo) {
      const url = `${appBase}/encyclopedia/#/material/${data.encyclopedia.material.material_id}`
      return (
        <Tooltip title="Show the material of this entry in the NOMAD Encyclopedia.">
          <Link href={url}>{name}</Link>
        </Tooltip>
      )
    } else {
      return name
    }
  }

  if (!exists) {
    return <EntryPageContent className={classes.root} fixed>
      <Typography className={classes.error}>
        This entry does not exist.
      </Typography>
    </EntryPageContent>
  }

  return (
    <Box className={classes.root} padding={3} margin="auto">
      <Grid container spacing={2}>

        <Grid item xs={4}>
          <Card className={classes.topCard}>
            <CardHeader
              title="Material"
              action={<ApiDialogButton title="Repository JSON" data={calcData} />}
            />
            <CardContent classes={{root: classes.cardContent}}>
              <Quantity column>
                <Quantity row>
                  <Quantity quantity="formula" label='formula' noWrap data={repo}/>
                  <Quantity quantity="dft.system" label='material type' noWrap data={repo}/>
                  <Quantity quantity={material_name} label='material name' noWrap data={repo}/>
                </Quantity>
                <Quantity row>
                  <Quantity quantity="dft.crystal_system" label='crystal system' noWrap data={repo}/>
                  <Quantity quantity="dft.spacegroup_symbol" label="spacegroup" noWrap data={repo}>
                    <Typography noWrap>
                      {normalizeDisplayValue(_.get(repo, 'dft.spacegroup_symbol'))} ({normalizeDisplayValue(_.get(repo, 'dft.spacegroup'))})
                    </Typography>
                  </Quantity>
                </Quantity>
              </Quantity>
              {reprSys
                ? <Structure
                  system={reprSys}
                ></Structure>
                : <Skeleton variant="rect" width={210} height={118} />
              }
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={4}>
          <Card className={classes.topCard}>
            <CardHeader
              title="Method"
              action={<ApiDialogButton title="Repository JSON" data={calcData} />}
            />
            <CardContent classes={{root: classes.cardContent}}>
              <Quantity row>
                <Quantity quantity="dft.code_name" label='dft code' noWrap data={repo}/>
                <Quantity quantity="dft.code_version" label='dft code version' noWrap data={repo}/>
              </Quantity>
              <Quantity row>
                <Quantity quantity="dft.basis_set" label='basis set' noWrap data={repo}/>
                <Quantity quantity="dft.xc_functional" label='xc functional' noWrap data={repo}/>
              </Quantity>
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={4}>
          <Card className={classes.topCard}>
            <CardHeader title="Entry" />
            <CardContent>
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
                <Quantity quantity="mainfile" loading={loadingRepo} noWrap ellipsisFront {...quantityProps} withClipboard />
                <Quantity quantity="calc_id" label={`${domain ? domain.entryLabel : 'entry'} id`} noWrap withClipboard {...quantityProps} />
                <Quantity quantity="upload_id" label='upload id' {...quantityProps} noWrap withClipboard />
                <Quantity quantity="upload_time" label='upload time' noWrap {...quantityProps} >
                  <Typography noWrap>
                    {new Date(calcData.upload_time).toLocaleString()}
                  </Typography>
                </Quantity>
                <Quantity quantity="raw_id" label='raw id' loading={loadingRepo} noWrap hideIfUnavailable {...quantityProps} withClipboard />
                <Quantity quantity="external_id" label='external id' loading={loadingRepo} hideIfUnavailable noWrap {...quantityProps} withClipboard />
              </Quantity>
            </CardContent>
          </Card>
        </Grid>

        {electronicStructure
          ? <Grid item xs={!!electronicStructure.bs * 7 + !!electronicStructure.dos * 5}>
            <Card>
              <CardHeader
                title="Electronic properties"
              />
              <CardContent classes={{root: classes.cardContent}}>
                <Box style={{width: '100%', height: '35rem', margin: '0 1rem'}}>
                  <ElectronicStructureOverview
                    data={electronicStructure}>
                  </ElectronicStructureOverview>
                </Box>
              </CardContent>
            </Card>
          </Grid>
          : null
        }
      </Grid>
    </Box>
  )
}

OverviewView.propTypes = {
  uploadId: PropTypes.string.isRequired,
  calcId: PropTypes.string.isRequired
}
