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
import { useRecoilValue } from 'recoil'
import { Box, Card, CardContent, Grid, Typography, Link, makeStyles, Divider } from '@material-ui/core'
import _ from 'lodash'
import { apiContext as apiContextV0 } from '../api'
import { useApi } from '../apiV1'
import { ApiDialog } from '../ApiDialogButton'
import ElectronicProperties from '../visualization/ElectronicProperties'
import VibrationalProperties from '../visualization/VibrationalProperties'
import GeometryOptimization from '../visualization/GeometryOptimization'
import Structure from '../visualization/Structure'
import Actions from '../Actions'
import Quantity from '../Quantity'
import { Link as RouterLink } from 'react-router-dom'
import { DOI } from '../search/DatasetList'
import { errorContext } from '../errors'
import {
  authorList,
  convertSI,
  getHighestOccupiedEnergy,
  toMateriaStructure,
  mergeObjects
} from '../../utils'
import { unitsState } from '../archive/ArchiveBrowser'
import { resolveRef, refPath } from '../archive/metainfo'
import searchQuantities from '../../searchQuantities'

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
  },
  content: {
    paddingBottom: 16,
    '&:last-child': {
      paddingBottom: 16
    }
  }
}))

function PropertyCard({title, children, actions}) {
  const styles = usePropertyCardStyles()
  return <Card className={styles.root}>
    <CardContent classes={{root: styles.content}}>
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
  const styles = useSidebarCardStyles()
  return <CardContent className={styles.content}>
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
const DFTEntryOverview = ({data}) => {
  // Determine which information source will be used: section_results of
  // section_metadata
  const hasResults = !!data?.results
  const units = useRecoilValue(unitsState)

  // Determine the set of available properties.
  const availableProps = useMemo(() => {
    let properties
    if (hasResults) {
      if (data?.results?.properties?.available_properties) {
        properties = new Set(data?.results?.properties?.available_properties)
      } else {
        properties = new Set()
      }
    } else {
      properties = new Set()
      properties.add('structure_original')
      if (data?.dft?.searchable_quantities) {
        const oldProps = new Set(data.dft.searchable_quantities)
        if (oldProps.has('electronic_dos')) {
          properties.add('dos_electronic')
        }
        if (oldProps.has('electronic_band_structure')) {
          properties.add('band_structure_electronic')
        }
        if (oldProps.has('phonon_dos')) {
          properties.add('dos_phonon')
        }
        if (oldProps.has('phonon_band_structure')) {
          properties.add('band_structure_phonon')
        }
        if (oldProps.has('thermodynamical_property_heat_capacity_C_v')) {
          properties.add('heat_capacity_constant_volume')
        }
        if (oldProps.has('vibrational_free_energy_at_constant_volume')) {
          properties.add('energy_free_helmholtz')
        }
        if (oldProps.has('vibrational_free_energy_at_constant_volume')) {
          properties.add('energy_free_helmholtz')
        }
      }
      if (data?.dft?.workflow?.workflow_type === 'geometry_optimization') {
        properties.add('geometry_optimization')
      }
    }
    return properties
  }, [data, hasResults])

  const apiV0 = useContext(apiContextV0).api
  const apiV1 = useApi()
  const {raiseError} = useContext(errorContext)
  const [dosElectronic, setDosElectronic] = useState(availableProps.has('dos_electronic') ? null : false)
  const [bsElectronic, setBsElectronic] = useState(availableProps.has('band_structure_electronic') ? null : false)
  const [bsPhonon, setBsPhonon] = useState(availableProps.has('band_structure_phonon') ? null : false)
  const [dosPhonon, setDosPhonon] = useState(availableProps.has('dos_phonon') ? null : false)
  const [heatCapacity, setHeatCapacity] = useState(availableProps.has('heat_capacity_constant_volume') ? null : false)
  const [freeEnergy, setFreeEnergy] = useState(availableProps.has('energy_free_helmholtz') ? null : false)
  const [geoOpt, setGeoOpt] = useState(availableProps.has('geometry_optimization') ? null : false)
  const [structures, setStructures] = useState(availableProps.has('structure_original') ? null : false)
  const [method, setMethod] = useState(null)
  const [loading, setLoading] = useState(true)
  const [showAPIDialog, setShowAPIDialog] = useState(false)
  const styles = useStyles()

  // Load the archive data from the server. If section_results is present, only
  // part of the archive will be loaded. Otherwise a full archive download is
  // performed.
  useEffect(() => {
    if (hasResults) {
      apiV1.results(data.entry_id).then(archive => {
        const url = `/entry/id/${data.upload_id}/${data.calc_id}/archive`
        const urlStruct = `${url}/results/properties/structures`

        // Structures
        let structs = []
        const original = toMateriaStructure(
          archive?.results?.properties?.structures?.structure_original,
          'original',
          `${urlStruct}/structure_original`
        )
        original && structs.push(original)
        const conventional = toMateriaStructure(
          archive?.results?.properties?.structures?.structure_conventional,
          'conventional',
          `${urlStruct}/structure_conventional`
        )
        conventional && structs.push(conventional)
        const primitive = toMateriaStructure(
          archive?.results?.properties?.structures?.structure_primitive,
          'primitive',
          `${urlStruct}/structure_primitive`
        )
        primitive && structs.push(primitive)
        structs = structs.length > 0 ? structs : false
        setStructures(structs)

        // Electronic properties
        const electronicDOS = archive?.results?.properties?.electronic?.dos_electronic
        if (electronicDOS) {
          const channel_info = electronicDOS.channel_info
          setDosElectronic({
            energies: resolveRef(electronicDOS.energies, archive),
            densities: resolveRef(electronicDOS.densities, archive),
            energy_highest_occupied: channel_info && Math.max(
              ...channel_info.map(x => x.energy_highest_occupied)
            ),
            m_path: `${url}/${refPath(electronicDOS.energies.split('/').slice(0, -1).join('/'))}`
          })
        }
        const electronicBS = archive?.results?.properties?.electronic?.band_structure_electronic
        if (electronicBS) {
          const channel_info = electronicBS.channel_info
          setBsElectronic({
            reciprocal_cell: resolveRef(electronicBS.reciprocal_cell, archive),
            segments: resolveRef(electronicBS.segments, archive),
            energy_highest_occupied: channel_info && Math.max(
              ...channel_info.map(x => x.energy_highest_occupied)
            ),
            channel_info: channel_info,
            m_path: `${url}/${refPath(electronicBS.reciprocal_cell.split('/').slice(0, -1).join('/'))}`
          })
        }
        // Geometry optimization
        const geoOptProps = archive?.results?.properties?.geometry_optimization
        const geoOptMethod = archive?.results?.method?.geometry_optimization
        if (geoOptProps) {
          setGeoOpt({
            energies: resolveRef(geoOptProps.energies, archive),
            energy_change_criteria: geoOptMethod?.input_energy_difference_tolerance
          })
        }
        // Vibrational properties
        const dosPhononProp = archive?.results?.properties?.vibrational?.dos_phonon
        const bsPhononProp = archive?.results?.properties?.vibrational?.band_structure_phonon
        const energyFreeProp = archive?.results?.properties?.vibrational?.energy_free_helmholtz
        const heatCapacityProp = archive?.results?.properties?.vibrational?.heat_capacity_constant_volume
        dosPhononProp && setDosPhonon({
          energies: resolveRef(dosPhononProp.energies, archive),
          densities: resolveRef(dosPhononProp.densities, archive),
          m_path: `${url}/${refPath(dosPhononProp.energies.split('/').slice(0, -1).join('/'))}`
        })
        bsPhononProp && setBsPhonon({
          segments: resolveRef(bsPhononProp.segments, archive),
          m_path: `${url}/${refPath(bsPhononProp.segments[0].split('/').slice(0, -1).join('/'))}`
        })
        energyFreeProp && setFreeEnergy({
          energies: resolveRef(energyFreeProp.energies, archive),
          temperatures: resolveRef(energyFreeProp.temperatures, archive),
          m_path: `${url}/${refPath(energyFreeProp.temperatures.split('/').slice(0, -1).join('/'))}`
        })
        heatCapacityProp && setHeatCapacity({
          heat_capacities: resolveRef(heatCapacityProp.heat_capacities, archive),
          temperatures: resolveRef(heatCapacityProp.temperatures, archive),
          m_path: `${url}/${refPath(heatCapacityProp.temperatures.split('/').slice(0, -1).join('/'))}`
        })
      }).catch(error => {
        if (error.name === 'DoesNotExist') {
        } else {
          raiseError(error)
        }
      })
    } else {
      apiV0.archive(data.upload_id, data.calc_id).then(archive => {
        let structs = []
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
            const doses = scc.section_dos
            const bss = scc.section_k_band
            if (!e_dos && doses) {
              for (let j = doses.length - 1; j > -1; --j) {
                const dos = doses[j]
                if (dos.dos_kind !== 'vibrational') {
                  e_dos = {
                    energies: dos.dos_energies_normalized,
                    densities: dos.dos_values_normalized,
                    energy_highest_occupied: 0,
                    m_path: `${url}/section_run/section_single_configuration_calculation:${i}/section_dos:${j}`
                  }
                }
              }
            }
            if (!e_bs && bss) {
              for (let j = bss.length - 1; j > -1; --j) {
                const band = scc.section_k_band[j]
                if (band.band_structure_kind !== 'vibrational') {
                  let channel_info
                  if (band.section_band_gap) {
                    channel_info = []
                    for (let k = 0; k < band.section_band_gap.length; ++k) {
                      const gap = band.section_band_gap[k]
                      channel_info.push({
                        index: k,
                        band_gap: gap.value,
                        band_gap_type: gap.type
                      })
                    }
                  }
                  e_bs = {
                    segments: band.section_k_band_segment,
                    reciprocal_cell: band.reciprocal_cell,
                    energy_highest_occupied: getHighestOccupiedEnergy(band, scc),
                    channel_info: channel_info,
                    m_path: `${url}/section_run/section_single_configuration_calculation:${i}/section_k_band:${j}`
                  }
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
                energies.push(e)
              }
            } catch (err) {
              failed = true
            }
            if (!failed) {
              const e_criteria_wf = section_wf?.section_geometry_optimization?.input_energy_difference_tolerance
              const sampling_method = section_run?.section_sampling_method
              const e_criteria_fs = sampling_method && sampling_method[0]?.geometry_optimization_energy_change
              const e_criteria = convertSI(e_criteria_wf || e_criteria_fs, 'joule', {energy: 'electron_volt'}, false)
              setGeoOpt({energies: energies, energy_change_criteria: e_criteria})
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
                segments: scc.section_k_band[scc.section_k_band.length - 1].section_k_band_segment,
                m_path: `${url}/${refPath(scc_ref)}/section_k_band:${scc.section_k_band.length - 1}`
              }
              v_dos = {
                energies: scc.section_dos[scc.section_dos.length - 1].dos_energies,
                densities: scc.section_dos[scc.section_dos.length - 1].dos_values_normalized,
                m_path: `${url}/${refPath(scc_ref)}/section_dos:${scc.section_dos.length - 1}`
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
                  heat_capacities: properties.thermodynamical_property_heat_capacity_C_v,
                  temperatures: properties.thermodynamical_property_temperature,
                  m_path: `${url}/section_run/section_frame_sequence:${sequences.length - 1}/section_thermodynamical_properties/thermodynamical_property_heat_capacity_C_v`
                }
                free_energy = {
                  energies: properties.vibrational_free_energy_at_constant_volume,
                  temperatures: properties.thermodynamical_property_temperature,
                  m_path: `${url}/section_run/section_frame_sequence:${sequences.length - 1}/section_thermodynamical_properties/vibrational_free_energy_at_constant_volume`
                }
              }
            }
            v_dos && setDosPhonon(v_dos)
            v_bs && setBsPhonon(v_bs)
            free_energy && setFreeEnergy(free_energy)
            heat_capacity && setHeatCapacity(heat_capacity)
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
            method_name: es_method,
            van_der_Waals_method: vdw_method,
            relativity_method: relativity_method,
            basis_set_name: basis_set
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
                m_path: `${url}/section_run/section_system:${i}`,
                name: 'original'
              }
              structs.push(reprSys)
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
            m_path: `${url}/section_metadata/encyclopedia/material/idealized_structure`,
            name: 'conventional'
          }
          structs.push(ideal)
        }
        setStructures(structs)
      }).catch(error => {
        if (error.name === 'DoesNotExist') {
        } else {
          raiseError(error)
        }
      }).finally(() => setLoading(false))
    }
  }, [data, apiV0, apiV1, raiseError, hasResults])

  const quantityProps = {data: data, loading: !data}
  const materialId = data?.results?.material?.material_id || data?.encyclopedia?.material?.material_id

  return (
    <Grid container spacing={0} className={styles.root}>

      {/* Left column */}
      <Grid item xs={4} className={styles.leftSidebar}>
        <SidebarCard title='Method'>
          <Quantity flex>
            {hasResults
              ? <>
                <Quantity quantity="results.method.simulation.program_name" label='program name' noWrap {...quantityProps}/>
                <Quantity quantity="results.method.simulation.program_version" label='program version' ellipsisFront {...quantityProps}/>
                <Quantity quantity="results.method.method_name" label='method name' noWrap {...quantityProps}/>
                {data?.results?.method?.method_name === 'DFT' && <>
                  <Quantity quantity="results.method.simulation.dft.xc_functional_type" label='xc functional family' noWrap {...quantityProps}/>
                  <Quantity quantity="results.method.simulation.dft.xc_functional_names" label='xc functional names' noWrap {...quantityProps}/>
                  <Quantity quantity="results.method.simulation.dft.basis_set_type" label='basis set type' noWrap {...quantityProps}/>
                  <Quantity quantity="results.method.simulation.dft.basis_set_name" label='basis set name' noWrap hideIfUnavailable {...quantityProps}/>
                  <Quantity quantity="results.method.simulation.dft.van_der_Waals_method" description="The used Van der Waals method." label='van der Waals method' noWrap hideIfUnavailable {...quantityProps}/>
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
              </>
              : <>
                <Quantity
                  quantity="dft.code_name"
                  description={searchQuantities['results.method.simulation.program_name']?.description || ''}
                  label='program name'
                  noWrap
                  {...quantityProps}
                />
                <Quantity
                  quantity="dft.code_version"
                  description={searchQuantities['results.method.simulation.program_version']?.description || ''}
                  label='program version'
                  ellipsisFront
                  {...quantityProps}
                />
                <Quantity
                  quantity="method_name"
                  description={searchQuantities['results.method.method_name']?.description || ''}
                  label='method name'
                  loading={loading}
                  noWrap
                  data={method}
                />
                <Quantity
                  quantity="dft.xc_functional"
                  description={searchQuantities['results.method.simulation.dft.xc_functional_type']?.description || ''}
                  label='xc functional family'
                  noWrap
                  {...quantityProps}
                />
                <Quantity
                  quantity="dft.xc_functional_names"
                  description={searchQuantities['results.method.simulation.dft.xc_functional_names']?.description || ''}
                  label='xc functional names'
                  noWrap
                  {...quantityProps}
                />
                <Quantity
                  quantity="dft.basis_set"
                  description={searchQuantities['results.method.simulation.dft.basis_set_type']?.description || ''}
                  label='basis set type'
                  noWrap
                  {...quantityProps}
                />
                <Quantity
                  quantity="basis_set_name"
                  description={searchQuantities['results.method.simulation.dft.basis_set_name']?.description || ''}
                  label='basis set name'
                  hideIfUnavailable
                  loading={loading}
                  noWrap
                  data={method}
                />
                <Quantity
                  quantity="van_der_Waals_method"
                  description="The used Van der Waals method."
                  label='van der Waals method'
                  hideIfUnavailable
                  loading={loading}
                  noWrap
                  data={method}
                />
                <Quantity
                  quantity="relativity_method"
                  description={searchQuantities['results.method.simulation.dft.relativity_method']?.description || ''}
                  label='relativity method'
                  hideIfUnavailable
                  loading={loading}
                  noWrap
                  data={method}
                />
              </>
            }
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
              description={searchQuantities['datasets'] && searchQuantities['datasets'].description}
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
            <Quantity quantity="entry_id" label='entry id' noWrap withClipboard {...quantityProps}/>
            { hasResults
              ? <Quantity quantity="results.material.material_id" label='material id' noWrap withClipboard {...quantityProps}/>
              : <Quantity
                quantity="encyclopedia.material.material_id"
                description={searchQuantities['results.material.material_id']?.description || ''}
                label='material id'
                noWrap
                withClipboard
                {...quantityProps}
              />
            }
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
                  {hasResults
                    ? <>
                      <Quantity quantity="results.material.chemical_formula_hill" label='formula' noWrap {...quantityProps}/>
                      <Quantity quantity="results.material.type_structural" label='structural type' noWrap {...quantityProps}/>
                      <Quantity quantity="results.material.material_name" label='material name' noWrap {...quantityProps}/>
                      {data?.results?.material?.symmetry &&
                      <Quantity row>
                        <Quantity
                          quantity="results.material.symmetry.crystal_system"
                          label='crystal system'
                          noWrap
                          {...quantityProps}
                        />
                        <Quantity
                          description="Space group symbol and number"
                          label="space group"
                          noWrap
                          {...quantityProps}
                        >
                          <Typography noWrap>
                            {normalizeDisplayValue(_.get(data, 'results.material.symmetry.space_group_symbol'))} ({normalizeDisplayValue(_.get(data, 'results.material.symmetry.space_group_number'))})
                          </Typography>
                        </Quantity>
                      </Quantity>
                      }
                    </>
                    : <>
                      <Quantity
                        quantity="formula"
                        description={searchQuantities['results.material.chemical_formula_hill']?.description || ''}
                        label='formula'
                        noWrap
                        {...quantityProps}
                      />
                      <Quantity
                        quantity="dft.system"
                        description={searchQuantities['results.material.type_structural']?.description || ''}
                        label='structural type'
                        noWrap
                        {...quantityProps}
                      />
                      <Quantity
                        quantity="encyclopedia.material.material_name"
                        description={searchQuantities['results.material.material_name']?.description || ''}
                        label='material name'
                        noWrap
                        {...quantityProps}
                      />
                      {data?.encyclopedia?.material?.bulk &&
                      <Quantity row>
                        <Quantity
                          quantity="dft.crystal_system"
                          description={searchQuantities['results.material.symmetry.crystal_system']?.description || ''}
                          label='crystal system'
                          noWrap
                          {...quantityProps}
                        />
                        <Quantity
                          description="Space group symbol and number"
                          label="space group"
                          noWrap
                          {...quantityProps}
                        >
                          <Typography noWrap>
                            {normalizeDisplayValue(_.get(data, 'dft.spacegroup_symbol'))} ({normalizeDisplayValue(_.get(data, 'dft.spacegroup'))})
                          </Typography>
                        </Quantity>
                      </Quantity>
                      }
                    </>
                  }
                </Quantity>
                {encyclopediaEnabled && materialId &&
                  <Actions
                    className={styles.actions}
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
                materialType={data?.results?.material?.type_structural || data?.dft?.system}
                aspectRatio={1.5}
                data-testid="viewer-material"
              />
            </Grid>
          </Grid>
        </PropertyCard>
        {(dosElectronic !== false ||
          bsElectronic !== false) &&
          <PropertyCard title="Electronic properties">
            <ElectronicProperties
              bs={bsElectronic}
              dos={dosElectronic}
              units={units}
            />
          </PropertyCard>
        }
        {geoOpt !== false &&
          <PropertyCard title="Geometry optimization">
            <GeometryOptimization data={geoOpt} units={units}/>
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
              units={units}
            />
          </PropertyCard>
        }
      </Grid>
    </Grid>
  )
}

DFTEntryOverview.propTypes = {
  data: PropTypes.object.isRequired
}

export default DFTEntryOverview
