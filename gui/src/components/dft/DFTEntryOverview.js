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
import { Box, Card, CardContent, Grid, Typography, Link, makeStyles, Divider } from '@material-ui/core'
import { apiContext } from '../api'
import ElectronicStructureOverview from '../visualization/ElectronicStructureOverview'
import VibrationalOverview from '../visualization/VibrationalOverview'
import { ApiDialog } from '../ApiDialogButton'
import Structure from '../visualization/Structure'
import NoData from '../visualization/NoData'
import Actions from '../Actions'
import Quantity from '../Quantity'
import { RecoilRoot } from 'recoil'
import { Link as RouterLink } from 'react-router-dom'
import { DOI } from '../search/DatasetList'
import { domains } from '../domains'
import { errorContext } from '../errors'
import { authorList, convertSI, mergeObjects } from '../../utils'
import { resolveRef, refPath } from '../archive/metainfo'
import _ from 'lodash'

import {appBase, encyclopediaEnabled, normalizeDisplayValue} from '../../config'
import GeoOptOverview from '../visualization/GeoOptOverview'

const useHeaderStyles = makeStyles(theme => ({
  root: {
    marginBottom: theme.spacing(1),
    display: 'flex',
    flexDirection: 'row',
    alignContent: 'flex-start'
  },
  title: {
    fontSize: '1.25rem'
  }
}))

function Header({title, actions}) {
  const styles = useHeaderStyles()
  return <Box className={styles.root}>
    <Box flexGrow={1}>
      <Typography className={styles.title}>
        {title}
      </Typography>
    </Box>
    {actions}
  </Box>
}

Header.propTypes = {
  title: PropTypes.string,
  actions: PropTypes.any
}

const usePropertyCardStyles = makeStyles(theme => ({
  root: {
    marginBottom: theme.spacing(2)
  },
  title: {
    fontSize: '1.25rem'
  }
}))

function PropertyCard({title, children, actions}) {
  const classes = usePropertyCardStyles()
  return <Card className={classes.root}>
    <CardContent>
      <Header title={title} actions={actions}></Header>
      {children}
    </CardContent>
  </Card>
}

PropertyCard.propTypes = {
  children: PropTypes.any,
  actions: PropTypes.any,
  title: PropTypes.string
}

const useSidebarCardStyles = makeStyles(theme => ({
  content: {
    padding: 0,
    paddingBottom: theme.spacing(2)
  },
  title: {
    fontSize: '1.25rem',
    marginBottom: theme.spacing(1)
  }
}))

function SidebarCard({title, actions, children}) {
  const classes = useSidebarCardStyles()
  return <CardContent className={classes.content}>
    {(title || actions) && <Header title={title} actions={actions}></Header>}
    {children}
  </CardContent>
}

SidebarCard.propTypes = {
  title: PropTypes.string,
  actions: PropTypes.any,
  children: PropTypes.any
}

const useStyles = makeStyles(theme => ({
  root: {
    marginBottom: theme.spacing(4)
  },
  error: {
    marginTop: theme.spacing(2)
  },
  cardHeader: {
    paddingBottom: 0
  },
  sidebar: {
    paddingRight: theme.spacing(3)
  },
  divider: {
    marginTop: theme.spacing(1),
    marginBottom: theme.spacing(1)
  },
  materialText: {
    height: '100%',
    display: 'flex',
    justifyContent: 'space-between',
    flexDirection: 'column'
  }
}))

/**
 * Shows an informative overview about the selected entry.
 */
export default function DFTEntryOverview({data}) {
  const {api} = useContext(apiContext)
  const {raiseError} = useContext(errorContext)
  const [electronicStructure, setElectronicStructure] = useState(null)
  const [vibrationalData, setVibrationalData] = useState(null)
  const [geoOpt, setGeoOpt] = useState(null)
  const [structures, setStructures] = useState(null)
  const [method, setMethod] = useState(null)
  const [loading, setLoading] = useState(true)
  const [showAPIDialog, setShowAPIDialog] = useState(false)
  const materialType = data?.encyclopedia?.material?.material_type

  const styles = useStyles()

  // When loaded for the first time, start downloading the archive. Once
  // finished, determine the final layout based on it's contents. TODO: When we
  // have more information stored in the ES index, it can be used to select
  // which parts of the Archive should be downloaded to reduce bandwidth.
  useEffect(() => {
    api.archive(data.upload_id, data.calc_id).then(archive => {
      let structs = {}
      const url = `/entry/id/${data.upload_id}/${data.calc_id}/archive`

      // Figure out what properties are present by looping over the SCCS. This
      // information will eventually be directly available in the ES index.
      let e_dos = null
      let e_bs = null
      const section_run = archive.section_run[0]
      let section_method = null
      const sccs = section_run.section_single_configuration_calculation
      if (sccs) {
        for (let i = sccs.length - 1; i > -1; --i) {
          const scc = sccs[i]
          if (!e_dos && scc.section_dos) {
            const first_dos = scc.section_dos[scc.section_dos.length - 1]
            if (first_dos.dos_kind !== 'vibrational') {
              e_dos = {
                'section_system': scc.single_configuration_calculation_to_system_ref,
                'section_method': scc.single_configuration_calculation_to_system_ref,
                'section_dos': scc.section_dos[scc.section_dos.length - 1],
                'path': `${url}/section_run/section_single_configuration_calculation:${i}/section_dos:${scc.section_dos.length - 1}`
              }
            }
          }
          if (!e_bs && scc.section_k_band) {
            const first_band = scc.section_k_band[scc.section_k_band.length - 1]
            if (first_band.band_structure_kind !== 'vibrational') {
              e_bs = {
                'section_system': scc.single_configuration_calculation_to_system_ref,
                'section_method': scc.single_configuration_calculation_to_system_ref,
                'section_k_band': scc.section_k_band[scc.section_k_band.length - 1],
                'path': `${url}/section_run/section_single_configuration_calculation:${i}/section_k_band:${scc.section_k_band.length - 1}`
              }
            }
          }
          if (section_method !== false) {
            let iMethod = scc.single_configuration_to_calculation_method_ref
            if (section_method === null) {
              section_method = iMethod
            } else if (iMethod !== section_method) {
              section_method = false
            }
          }
        }
      }

      if (e_dos || e_bs) {
        setElectronicStructure({
          'dos': e_dos, 'bs': e_bs
        })
      }

      // See if there are workflow results
      const section_wf = archive.section_workflow
      if (section_wf) {
        const wfType = section_wf.workflow_type

        // Gather energies, trajectory and energy change threshold from geometry
        // optimization
        if (wfType === 'geometry_optimization') {
          let failed = false
          let energies = []
          const trajectory = []
          try {
            const calculations = section_wf.calculations_ref
            let initialEnergy = null
            if (!calculations) {
              throw Error('no calculations')
            }
            for (let i = 0; i < calculations.length; ++i) {
              let ref = calculations[i]
              const calc = resolveRef(ref, archive)
              let e = calc.energy_total
              if (e === undefined) {
                if (i === calculations.length - 1) {
                  break
                } else {
                  throw Error('invalid energy value')
                }
              }
              if (i === 0) {
                initialEnergy = e
              }
              energies.push(e - initialEnergy)
              let sys = calc.single_configuration_calculation_to_system_ref
              sys = resolveRef(sys, archive)
              if (sys === undefined) {
                throw Error('invalid system reference')
              }
              trajectory.push({
                species: sys.atom_species,
                cell: sys.lattice_vectors ? convertSI(sys.lattice_vectors, 'meter', {length: 'angstrom'}, false) : undefined,
                positions: convertSI(sys.atom_positions, 'meter', {length: 'angstrom'}, false),
                pbc: sys.configuration_periodic_dimensions
              })
            }
          } catch (err) {
            failed = true
          }
          if (!failed) {
            energies = convertSI(energies, 'joule', {energy: 'electron_volt'}, false)
            const e_criteria_wf = section_wf?.section_geometry_optimization?.input_energy_difference_tolerance
            const sampling_method = section_run?.section_sampling_method
            const e_criteria_fs = sampling_method && sampling_method[0]?.geometry_optimization_energy_change
            const e_criteria = e_criteria_wf || e_criteria_fs
            setGeoOpt({energies: energies, structures: trajectory, energy_change_criteria: e_criteria})
          } else {
            setGeoOpt({})
          }
        } else if (wfType === 'phonon') {
          // Find phonon dos and dispersion
          const scc_ref = section_wf.calculation_result_ref
          const scc = resolveRef(scc_ref, archive)
          let v_dos = null
          let v_bs = null
          if (scc) {
            v_bs = {
              section_system: scc.single_configuration_calculation_to_system_ref,
              section_method: scc.single_configuration_calculation_to_system_ref,
              section_k_band: scc.section_k_band[scc.section_k_band.length - 1],
              path: `${url}/${refPath(scc_ref)}/section_k_band:${scc.section_k_band.length - 1}`
            }
            v_dos = {
              section_system: scc.single_configuration_calculation_to_system_ref,
              section_method: scc.single_configuration_calculation_to_system_ref,
              section_dos: scc.section_dos[scc.section_dos.length - 1],
              path: `${url}/${refPath(scc_ref)}/section_dos:${scc.section_dos.length - 1}`
            }
          }

          // Find thermal properties
          let free_energy = null
          let heat_capacity = null
          let temperature = null
          const sequences = section_run.section_frame_sequence
          const sequence = sequences && sequences[sequences.length - 1]
          if (sequence) {
            const properties = sequence.section_thermodynamical_properties && sequence.section_thermodynamical_properties[0]
            if (properties) {
              heat_capacity = {
                thermodynamical_property_heat_capacity_C_v: properties.thermodynamical_property_heat_capacity_C_v,
                path: `${url}/section_run/section_frame_sequence:${sequences.length - 1}/section_thermodynamical_properties/thermodynamical_property_heat_capacity_C_v`
              }
              free_energy = {
                vibrational_free_energy_at_constant_volume: properties.vibrational_free_energy_at_constant_volume,
                path: `${url}/section_run/section_frame_sequence:${sequences.length - 1}/section_thermodynamical_properties/vibrational_free_energy_at_constant_volume`
              }
              temperature = properties.thermodynamical_property_temperature
            }
          }

          if (v_dos || v_bs || free_energy || heat_capacity) {
            setVibrationalData({
              dos: v_dos,
              bs: v_bs,
              free_energy: free_energy,
              heat_capacity: heat_capacity,
              temperature: temperature
            })
          }
        }
      }

      // Get method details. Any referenced core_setttings will also be taken
      // into account. If there were no SCCs from which the method could be
      // selected, simply select the last available method.
      if (section_method) {
        section_method = resolveRef(section_method, archive)
      } else {
        const methods = section_run.section_method
        if (methods) {
          section_method = methods[methods.length - 1]
        }
      }
      if (section_method) {
        const refs = section_method?.section_method_to_method_refs
        if (refs) {
          for (const ref of refs) {
            if (ref.method_to_method_kind === 'core_settings') {
              section_method = mergeObjects(resolveRef(ref.method_to_method_ref, archive), section_method)
            }
          }
        }
        const es_method = section_method?.electronic_structure_method
        const vdw_method = section_method?.van_der_Waals_method
        const relativity_method = section_method?.relativity_method
        const basis_set = section_method?.basis_set
        setMethod({
          electronic_structure_method: es_method,
          van_der_Waals_method: vdw_method,
          relativity_method: relativity_method,
          basis_set: basis_set
        })
      }

      // Get the representative system by looping over systems
      let reprSys = null
      const systems = section_run.section_system
      if (systems) {
        for (let i = systems.length - 1; i > -1; --i) {
          const sys = systems[i]
          if (!reprSys && sys.is_representative) {
            const reprSys = {
              species: sys.atom_species,
              cell: sys.lattice_vectors ? convertSI(sys.lattice_vectors, 'meter', {length: 'angstrom'}, false) : undefined,
              positions: convertSI(sys.atom_positions, 'meter', {length: 'angstrom'}, false),
              pbc: sys.configuration_periodic_dimensions,
              path: `${url}/section_run/section_system:${i}`
            }
            structs.original = reprSys
            break
          }
        }
      }

      // Get the conventional (=normalized) system, if present
      let idealSys = archive?.section_metadata?.encyclopedia?.material?.idealized_structure
      if (idealSys && data?.dft?.system === 'bulk') {
        const ideal = {
          species: idealSys.atom_labels,
          cell: idealSys.lattice_vectors ? convertSI(idealSys.lattice_vectors, 'meter', {length: 'angstrom'}, false) : undefined,
          positions: idealSys.atom_positions,
          fractional: true,
          pbc: idealSys.periodicity,
          path: `${url}/section_metadata/encyclopedia/material/idealized_structure`
        }
        structs.conventional = ideal
      }
      setStructures(structs)
      setLoading(false)
    }).catch(error => {
      setLoading(false)
      if (error.name === 'DoesNotExist') {
      } else {
        raiseError(error)
      }
    })
  }, [data, api, raiseError, setElectronicStructure, setStructures])

  const quantityProps = {data: data, loading: !data}
  const domain = data.domain && domains[data.domain]

  return (
    <RecoilRoot>
      <Grid container spacing={0} className={styles.root}>

        {/* Left column */}
        <Grid item xs={4} className={styles.sidebar}>
          <SidebarCard title='Method'>
            <Quantity flex>
              <Quantity quantity="dft.code_name" label='code name' noWrap {...quantityProps}/>
              <Quantity quantity="dft.code_version" label='code version' noWrap {...quantityProps}/>
              <Quantity
                quantity="electronic_structure_method"
                label='electronic structure method'
                loading={loading}
                description="The used electronic structure method."
                noWrap
                data={method}
              />
              <Quantity quantity="dft.xc_functional" label='xc functional family' noWrap {...quantityProps}/>
              <Quantity quantity="dft.xc_functional_names" label='xc functional names' noWrap {...quantityProps}/>
              <Quantity quantity="dft.basis_set" label='basis set type' noWrap {...quantityProps}/>
              <Quantity quantity="basis_set" label='basis set name' noWrap hideIfUnavailable data={method}/>
              {method?.van_der_Waals_method && <Quantity quantity="van_der_Waals_method" label='van der Waals method' noWrap {...quantityProps}/>}
              {method?.relativity_method && <Quantity quantity="relativity_method" label='relativity method' noWrap data={method}/>}
            </Quantity>
          </SidebarCard>
          <Divider className={styles.divider} />
          <SidebarCard title='Author metadata'>
            <Quantity flex>
              <Quantity quantity='comment' placeholder='no comment' {...quantityProps} />
              <Quantity quantity='references' placeholder='no references' {...quantityProps}>
                {data.references &&
                <div style={{display: 'inline-grid'}}>
                  {data.references.map(ref => <Typography key={ref} noWrap>
                    <Link href={ref}>{ref}</Link>
                  </Typography>)}
                </div>}
              </Quantity>
              <Quantity quantity='authors' {...quantityProps}>
                <Typography>
                  {authorList(data || [])}
                </Typography>
              </Quantity>
              <Quantity quantity='datasets' placeholder='no datasets' {...quantityProps}>
                {data.datasets &&
                <div>
                  {data.datasets.map(ds => (
                    <Typography key={ds.dataset_id}>
                      <Link component={RouterLink} to={`/dataset/id/${ds.dataset_id}`}>{ds.name}</Link>
                      {ds.doi ? <span>&nbsp; (<DOI doi={ds.doi}/>)</span> : ''}
                    </Typography>))}
                </div>}
              </Quantity>
            </Quantity>
          </SidebarCard>
          <Divider className={styles.divider}/>
          <SidebarCard>
            <Quantity column style={{maxWidth: 350}}>
              <Quantity quantity="mainfile" noWrap ellipsisFront withClipboard {...quantityProps}/>
              <Quantity quantity="calc_id" label={`${domain ? domain.entryLabel : 'entry'} id`} noWrap withClipboard {...quantityProps}/>
              <Quantity quantity="encyclopedia.material.material_id" label='material id' noWrap withClipboard {...quantityProps}/>
              <Quantity quantity="upload_id" label='upload id' noWrap withClipboard {...quantityProps}/>
              <Quantity quantity="upload_time" label='upload time' noWrap {...quantityProps}>
                <Typography noWrap>
                  {new Date(data.upload_time).toLocaleString()}
                </Typography>
              </Quantity>
              <Quantity quantity="raw_id" label='raw id' noWrap hideIfUnavailable withClipboard {...quantityProps}/>
              <Quantity quantity="external_id" label='external id' hideIfUnavailable noWrap withClipboard {...quantityProps}/>
              <Quantity quantity="last_processing" label='last processing' placeholder="not processed" noWrap {...quantityProps}>
                <Typography noWrap>
                  {new Date(data.last_processing).toLocaleString()}
                </Typography>
              </Quantity>
              <Quantity quantity="last_processing" label='processing version' noWrap placeholder="not processed" {...quantityProps}>
                <Typography noWrap>
                  {data.nomad_version}/{data.nomad_commit}
                </Typography>
              </Quantity>
            </Quantity>
          </SidebarCard>
          <ApiDialog data={data} open={showAPIDialog} onClose={() => { setShowAPIDialog(false) }}></ApiDialog>
          <Actions
            justifyContent='flex-end'
            variant='outlined'
            color='primary'
            size='medium'
            actions={[{
              tooltip: 'Show the API access code',
              onClick: (event) => { setShowAPIDialog(!showAPIDialog) },
              content: 'API'
            }]}
          >
          </Actions>
        </Grid>

        {/* Right column */}
        <Grid item xs={8}>
          <PropertyCard title="Material">
            <Grid container spacing={1}>
              <Grid item xs={5}>
                <Box className={styles.materialText}>
                  <Quantity column>
                    <Quantity quantity="formula" label='formula' noWrap {...quantityProps}/>
                    <Quantity quantity="dft.system" label='material type' noWrap {...quantityProps}/>
                    <Quantity quantity="encyclopedia.material.material_name" label='material name' noWrap {...quantityProps}/>
                    <Quantity row>
                      {materialType === 'bulk' && <Quantity quantity="dft.crystal_system" label='crystal system' noWrap {...quantityProps}/>}
                      {materialType === 'bulk' && <Quantity quantity="dft.spacegroup_symbol" label="spacegroup" noWrap {...quantityProps}>
                        <Typography noWrap>
                          {normalizeDisplayValue(_.get(data, 'dft.spacegroup_symbol'))} ({normalizeDisplayValue(_.get(data, 'dft.spacegroup'))})
                        </Typography>
                      </Quantity>}
                    </Quantity>
                  </Quantity>
                  {encyclopediaEnabled && data?.encyclopedia?.material?.material_id &&
                    <Actions
                      className={styles.actions}
                      justifyContent='flex-start'
                      color='primary'
                      variant='text'
                      size='medium'
                      actions={[{
                        tooltip: 'View this material in the Encyclopedia',
                        content: 'Encyclopedia',
                        href: `${appBase}/encyclopedia/#/material/${data.encyclopedia.material.material_id}`
                      }]}
                    >
                    </Actions>
                  }
                </Box>
              </Grid>
              <Grid item xs={7} style={{marginTop: '-2rem'}}>
                {(loading || !_.isEmpty(structures))
                  ? <Structure systems={structures} materialType={data?.dft?.system} aspectRatio={1.5}/>
                  : <NoData aspectRatio={1.5}/>
                }
              </Grid>
            </Grid>
          </PropertyCard>
          {electronicStructure &&
            <PropertyCard title="Electronic properties">
              <ElectronicStructureOverview
                data={electronicStructure}>
              </ElectronicStructureOverview>
            </PropertyCard>
          }
          {geoOpt && structures &&
            <PropertyCard title="Geometry optimization">
              <GeoOptOverview data={geoOpt}></GeoOptOverview>
            </PropertyCard>
          }
          {vibrationalData &&
            <PropertyCard title="Vibrational properties">
              <VibrationalOverview
                data={vibrationalData}>
              </VibrationalOverview>
            </PropertyCard>
          }
        </Grid>
      </Grid>
    </RecoilRoot>
  )
}

DFTEntryOverview.propTypes = {
  data: PropTypes.object.isRequired
}
