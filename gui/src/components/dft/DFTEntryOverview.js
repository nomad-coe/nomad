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
import React, { useContext, useState, useMemo, useEffect } from 'react'
import PropTypes from 'prop-types'
import { Box, Card, CardContent, Grid, Typography, Link, makeStyles, Divider } from '@material-ui/core'
import { apiContext } from '../api'
import { apiContext as apiContextv1 } from '../apiv1'
import { ApiDialog } from '../ApiDialogButton'
import ElectronicProperties from '../visualization/ElectronicProperties'
import VibrationalProperties from '../visualization/VibrationalProperties'
import GeometryOptimization from '../visualization/GeometryOptimization'
import Structure from '../visualization/Structure'
import NoData from '../visualization/NoData'
import Actions from '../Actions'
import Quantity from '../Quantity'
import { RecoilRoot } from 'recoil'
import { Link as RouterLink } from 'react-router-dom'
import { DOI } from '../search/DatasetList'
import { errorContext } from '../errors'
import { authorList, convertSI, getHighestOccupiedEnergy } from '../../utils'
import { resolveRef, refPath } from '../archive/metainfo'
import searchQuantities from '../../searchQuantities'
import _ from 'lodash'

import {appBase, encyclopediaEnabled, normalizeDisplayValue} from '../../config'

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
  leftSidebar: {
    maxWidth: '32%',
    flexBasis: '32%',
    flexGrow: 0,
    paddingRight: theme.spacing(3)
  },
  rightSidebar: {
    maxWidth: '67.99%',
    flexBasis: '67.99%',
    flexGrow: 0
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
  const availableProps = useMemo(() => {
    let properties
    if (data?.dft?.searchable_quantities) {
      properties = new Set(data.dft.searchable_quantities)
    } else {
      properties = new Set()
    }
    if (data?.dft?.workflow?.workflow_type === 'geometry_optimization') {
      properties.add('geometry_optimization')
    }
    return properties
  }, [data])
  const {api} = useContext(apiContext)
  const apiv1 = useContext(apiContextv1).api
  const {raiseError} = useContext(errorContext)
  const [dosElectronic, setDosElectronic] = useState(availableProps.has('electronic_dos') ? null : false)
  const [bsElectronic, setBsElectronic] = useState(availableProps.has('electronic_band_structure') ? null : false)
  const [bsPhonon, setBsPhonon] = useState(availableProps.has('phonon_band_structure') ? null : false)
  const [dosPhonon, setDosPhonon] = useState(availableProps.has('phonon_dos') ? null : false)
  const [heatCapacity, setHeatCapacity] = useState(availableProps.has('thermodynamical_property_heat_capacity_C_v') ? null : false)
  const [freeEnergy, setFreeEnergy] = useState(availableProps.has('vibrational_free_energy_at_constant_volume') ? null : false)
  const [dataGeoOpt, setDataGeoOpt] = useState(availableProps.has('geometry_optimization') ? null : false)
  const [structures, setStructures] = useState(null)
  const [loading, setLoading] = useState(true)
  const [showAPIDialog, setShowAPIDialog] = useState(false)
  const styles = useStyles()

  // Load section_results from the server
  useEffect(() => {
    apiv1.results(data.entry_id).then(archive => {
    }).catch(error => {
      setLoading(false)
      if (error.name === 'DoesNotExist') {
      } else {
        raiseError(error)
      }
    })
  }, [data, apiv1, raiseError, setStructures])

  // When loaded for the first time, start downloading the archive. Once
  // finished, determine the final layout based on it's contents. TODO: When we
  // have more information stored in the ES index, it can be used to select
  // which parts of the Archive should be downloaded to reduce bandwidth.
  useEffect(() => {
    api.archive(data.upload_id, data.calc_id).then(archive => {
      let structs = {}
      const url = `/entry/id/${data.upload_id}/${data.calc_id}/archive`

      // Check that at least one section run is available. If not, break the execution.
      let section_run = archive?.section_run
      if (section_run) {
        const nRuns = section_run.length
        section_run = section_run[nRuns - 1]
      } else {
        return
      }

      // Figure out what properties are present by looping over the SCCS. This
      // information will eventually be directly available in the ES index.
      let e_dos = null
      let e_bs = null
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
            first_band.energy_highest_occupied = getHighestOccupiedEnergy(first_band, scc)
            if (first_band.band_structure_kind !== 'vibrational') {
              e_bs = {
                'section_system': scc.single_configuration_calculation_to_system_ref,
                'section_method': scc.single_configuration_calculation_to_system_ref,
                'section_k_band': first_band,
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
      if (e_dos) {
        setDosElectronic(e_dos)
      }
      if (e_bs) {
        setBsElectronic(e_bs)
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
            }
          } catch (err) {
            failed = true
          }
          if (!failed) {
            energies = convertSI(energies, 'joule', {energy: 'electron_volt'}, false)
            const e_criteria_wf = section_wf?.section_geometry_optimization?.input_energy_difference_tolerance
            const sampling_method = section_run?.section_sampling_method
            const e_criteria_fs = sampling_method && sampling_method[0]?.geometry_optimization_energy_change
            const e_criteria = convertSI(e_criteria_wf || e_criteria_fs, 'joule', {energy: 'electron_volt'}, false)
            setDataGeoOpt({energies: energies, energy_change_criteria: e_criteria})
          } else {
            setDataGeoOpt({})
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
          const sequences = section_run.section_frame_sequence
          const sequence = sequences && sequences[sequences.length - 1]
          if (sequence) {
            const properties = sequence.section_thermodynamical_properties && sequence.section_thermodynamical_properties[0]
            if (properties) {
              heat_capacity = {
                thermodynamical_property_heat_capacity_C_v: properties.thermodynamical_property_heat_capacity_C_v,
                path: `${url}/section_run/section_frame_sequence:${sequences.length - 1}/section_thermodynamical_properties/thermodynamical_property_heat_capacity_C_v`,
                temperature: properties.thermodynamical_property_temperature
              }
              free_energy = {
                vibrational_free_energy_at_constant_volume: properties.vibrational_free_energy_at_constant_volume,
                path: `${url}/section_run/section_frame_sequence:${sequences.length - 1}/section_thermodynamical_properties/vibrational_free_energy_at_constant_volume`,
                temperature: properties.thermodynamical_property_temperature
              }
            }
          }
          v_dos && setDosPhonon(v_dos)
          v_bs && setBsPhonon(v_bs)
          free_energy && setFreeEnergy(free_energy)
          heat_capacity && setHeatCapacity(heat_capacity)
        }
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
    }).catch(error => {
      if (error.name === 'DoesNotExist') {
      } else {
        raiseError(error)
      }
    }).finally(() => setLoading(false))
  }, [data, api, raiseError])

  const quantityProps = {data: data, loading: !data}

  return (
    <RecoilRoot>
      <Grid container spacing={0} className={styles.root}>

        {/* Left column */}
        <Grid item xs={4} className={styles.leftSidebar}>
          <SidebarCard title='Method'>
            <Quantity flex>
              <Quantity quantity="results.method.simulation.program_name" label='program name' noWrap {...quantityProps}/>
              <Quantity quantity="results.method.simulation.program_version" label='program version' ellipsisFront {...quantityProps}/>
              <Quantity quantity="results.method.method_name" label='electronic structure method' noWrap {...quantityProps}/>
              {data?.results?.method?.method_name === 'DFT' && <>
                <Quantity quantity="results.method.simulation.dft.xc_functional_type" label='xc functional family' noWrap {...quantityProps}/>
                <Quantity quantity="results.method.simulation.dft.xc_functional_names" label='xc functional names' noWrap {...quantityProps}/>
                <Quantity quantity="results.method.simulation.dft.basis_set_type" label='basis set type' noWrap {...quantityProps}/>
                <Quantity quantity="results.method.simulation.dft.basis_set_name" label='basis set name' noWrap hideIfUnavailable {...quantityProps}/>
                <Quantity quantity="results.method.simulation.dft.van_der_Waals_method" label='van der Waals method' noWrap hideIfUnavailable {...quantityProps}/>
                <Quantity quantity="results.method.simulation.dft.relativity_method" label='relativity method' noWrap hideIfUnavailable {...quantityProps}/>
                <Quantity quantity="results.method.simulation.dft.core_electron_treatment" label='core electron treatment' noWrap {...quantityProps}/>
                <Quantity quantity="results.method.simulation.dft.spin_polarized" label='spin-polarized' hideIfUnavailable noWrap {...quantityProps}/>
                <Quantity quantity="results.method.simulation.dft.smearing_type" label='smearing type' noWrap hideIfUnavailable {...quantityProps}/>
                <Quantity quantity="results.method.simulation.dft.smearing_width" label='smearing width' noWrap hideIfUnavailable {...quantityProps}/>
              </>}
              {data?.results?.method?.method_name === 'GW' && <>
                <Quantity quantity="results.method.simulation.gw.gw_type" label='gw type' noWrap {...quantityProps}/>
                <Quantity quantity="results.method.simulation.gw.starting_point" label='ground state xc functional' noWrap {...quantityProps}/>
              </>}
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
              <Quantity
                desciption={searchQuantities['datasets'] && searchQuantities['datasets'].description}
                label='datasets'
                placeholder='no datasets'
                {...quantityProps}
              >
                {(data.datasets && data.datasets.length !== 0) &&
                <div>
                  {data.datasets.map(ds => (
                    <Typography key={ds.dataset_id}>
                      <Link component={RouterLink} to={`/dataset/id/${ds.dataset_id}`}>{ds.name}</Link>
                      {ds.doi ? <span>&nbsp;<DOI style={{display: 'inline'}} parentheses doi={ds.doi}/></span> : ''}
                    </Typography>))}
                </div>}
              </Quantity>
            </Quantity>
          </SidebarCard>
          <Divider className={styles.divider}/>
          <SidebarCard>
            <Quantity column style={{maxWidth: 350}}>
              <Quantity quantity="mainfile" noWrap ellipsisFront withClipboard {...quantityProps}/>
              <Quantity quantity="entry_id" label='entry_id' noWrap withClipboard {...quantityProps}/>
              <Quantity quantity="results.material.material_id" label='material id' noWrap withClipboard {...quantityProps}/>
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
              <Quantity description="Version used in the last processing" label='processing version' noWrap placeholder="not processed" {...quantityProps}>
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
        <Grid item xs={8} className={styles.rightSidebar}>
          <PropertyCard title="Material">
            <Grid container spacing={1}>
              <Grid item xs={5}>
                <Box className={styles.materialText}>
                  <Quantity column>
                    <Quantity quantity="results.material.chemical_formula_hill" label='formula' noWrap {...quantityProps}/>
                    <Quantity quantity="results.material.type_structural" label='structural type' noWrap {...quantityProps}/>
                    <Quantity quantity="results.material.material_name" label='material name' noWrap {...quantityProps}/>
                    {data?.results?.material?.symmetry &&
                      <Quantity row>
                        <Quantity
                          quantity="results.material.symmetry.crystal_system"
                          label='crystal system'
                          hideIfUnavailable
                          noWrap
                          {...quantityProps}
                        />
                        <Quantity
                          description="Space group symbol and number"
                          label="spacegroup"
                          hideIfUnavailable
                          noWrap
                          {...quantityProps}
                        >
                          <Typography noWrap>
                            {normalizeDisplayValue(_.get(data, 'results.material.symmetry.space_group_symbol'))} ({normalizeDisplayValue(_.get(data, 'results.material.symmetry.space_group_number'))})
                          </Typography>
                        </Quantity>
                      </Quantity>
                    }
                  </Quantity>
                  {encyclopediaEnabled && data?.results?.material?.material_id &&
                    <Actions
                      className={styles.actions}
                      justifyContent='flex-start'
                      color='primary'
                      variant='text'
                      size='medium'
                      actions={[{
                        tooltip: 'View this material in the Encyclopedia',
                        content: 'Encyclopedia',
                        href: `${appBase}/encyclopedia/#/material/${data.results.material.material_id}`
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
          {(dosElectronic !== false ||
            bsElectronic !== false) &&
            <PropertyCard title="Electronic properties">
              <ElectronicProperties
                bs={bsElectronic}
                dos={dosElectronic}
              >
              </ElectronicProperties>
            </PropertyCard>
          }
          {dataGeoOpt !== false &&
            <PropertyCard title="Geometry optimization">
              <GeometryOptimization data={dataGeoOpt}></GeometryOptimization>
            </PropertyCard>
          }
          {(dosPhonon !== false ||
            bsPhonon !== false ||
            heatCapacity !== false ||
            freeEnergy !== false) &&
            <PropertyCard title="Vibrational properties">
              <VibrationalProperties
                bs={bsPhonon}
                dos={dosPhonon}
                heatCapacity={heatCapacity}
                freeEnergy={freeEnergy}
              >
              </VibrationalProperties>
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
