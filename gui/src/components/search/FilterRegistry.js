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

import React, { useEffect, useState, useMemo } from 'react'
import PropTypes from 'prop-types'
import { useGlobalMetainfo, getQuantities } from '../archive/metainfo'
import { setToArray, DType, getSuggestions, getOptions, getDatatype, glob, parseQuantityName } from '../../utils'
import { searchQuantities, schemaSeparator, dtypeSeparator, yamlSchemaPrefix } from '../../config'
import { Filter, getEnumOptions } from './Filter'
import elementData from '../../elementData'
import { Typography, Box } from '@material-ui/core'
import { useErrors } from '../errors'

// Containers for filter information
export const defaultFilterData = {} // Stores data for each registered filter
const dtypeMap = {
  [DType.Int]: 'int',
  [DType.Float]: 'float',
  [DType.Timestamp]: 'datetime',
  [DType.String]: 'str',
  [DType.Enum]: 'str',
  [DType.Boolean]: 'bool'
}

/**
 * This function is used to register a new filter within the SearchContext.
 * Filters are entities that can be searched through the filter panel and the
 * search bar, and can be encoded in the URL. Notice that a filter in this
 * context does not necessarily have to correspond to a quantity in the
 * metainfo, altought typically this is preferrable.
 *
 * Only registered filters may be searched for. The registration must happen
 * before any components use the filters. This is because:
 *  - The initial aggregation results must be fetched before any components
 *  using the filter values are rendered.
 *  - Several components need to know the list of available filters (e.g. the
 *  search bar and the search panel). If filters are only registered during
 *  component initialization, it may already be too late to update other
 *  components.
 *
 * @param {string} name Name of the filter.
 * @param {obj} config Data object containing options for the filter.
 */
function saveFilter(name, config, parent) {
  if (defaultFilterData[name]) {
    throw Error(`Trying to register filter "${name}"" multiple times.`)
  }
  const def = searchQuantities[name]
  const {path: quantity, schema} = parseQuantityName(name)
  const newConf = {
    ...(config || {}),
    quantity,
    schema,
    name: config?.name || def?.name || name,
    aggregatable: def?.aggregatable
  }
  const data = defaultFilterData[name] || new Filter(def, newConf, parent)
  defaultFilterData[name] = data
  return data
}

/**
 * Used to register a filter value for an individual metainfo quantity or
 * section.
 */
function registerFilter(name, config, subQuantities) {
  const parent = saveFilter(name, config)
  if (subQuantities) {
    for (const subConfig of subQuantities) {
      const subname = `${name}.${subConfig.name}`
      saveFilter(subname, subConfig, parent)
    }
  }
}

/**
 * Used to register a filter that is based on a subset of quantity values.
 */
function registerFilterOptions(name, target, label, description, options) {
  const keys = Object.keys(options)
  registerFilter(
    name,
    {
      aggs: {
        terms: {
          set: ((target) => (config) => ({quantity: target, include: keys, ...config}))(target),
          get: (() => (agg) => agg)(target)
        }
      },
      value: {
        set: (newQuery, oldQuery, value) => {
          const data = newQuery[target] || new Set()
          value.forEach((item) => { data.add(item) })
          newQuery[`${target}:all`] = data
        }
      },
      dtype: DType.Enum,
      multiple: true,
      repeats: true,
      exclusive: false,
      queryMode: 'all',
      options: options,
      label: label,
      description: description
    }
  )
}

// Presets for different kind of quantities
const termQuantity = {aggs: {terms: {size: 5}}}
const termQuantityLarge = {aggs: {terms: {size: 10}}}
const termQuantityBool = {
  aggs: {terms: {size: 2, include: ['false', 'true']}},
  options: {
    false: {label: 'false'},
    true: {label: 'true'}
  }
}
const termQuantityNonExclusive = {aggs: {terms: {size: 5}}, exclusive: false}
const termQuantityAll = {aggs: {terms: {size: 5}}, exclusive: false, multiple: true, queryMode: 'all'}
const termQuantityAllNonExclusive = {...termQuantityNonExclusive, queryMode: 'all'}
const noAggQuantity = {}
const nestedQuantity = {}
const noQueryQuantity = {multiple: false, global: true}
const numberHistogramQuantity = {multiple: false, exclusive: false}

// Filters that directly correspond to a metainfo value
registerFilter(
  'results.material.structural_type',
  {
    ...termQuantity,
    scale: 'log',
    label: "Dimensionality",
    options: getEnumOptions('results.material.structural_type', ['not processed', 'unavailable'])
  }
)

registerFilter('results.material.functional_type', termQuantityNonExclusive)
registerFilter('results.material.compound_type', termQuantityNonExclusive)
registerFilter('results.material.material_name', termQuantity)
registerFilter('results.material.chemical_formula_hill', {...termQuantity, placeholder: "E.g. H2O2, C2H5Br"})
registerFilter('results.material.chemical_formula_iupac', {...termQuantity, placeholder: "E.g. GaAs, SiC", label: 'Chemical formula IUPAC'})
registerFilter('results.material.chemical_formula_reduced', {...termQuantity, placeholder: "E.g. H2NaO, ClNa"})
registerFilter('results.material.chemical_formula_anonymous', {...termQuantity, placeholder: "E.g. A2B, A3B2C2"})
registerFilter('results.material.n_elements', {...numberHistogramQuantity})
registerFilter('results.material.symmetry.bravais_lattice', termQuantity)
registerFilter('results.material.symmetry.crystal_system', termQuantity)
registerFilter(
  'results.material.symmetry.structure_name',
  {
    ...termQuantity,
    options: getEnumOptions('results.material.symmetry.structure_name', ['not processed', 'cubic perovskite'])
  }
)
registerFilter('results.material.symmetry.strukturbericht_designation', termQuantity)
registerFilter('results.material.symmetry.space_group_symbol', {...termQuantity, placeholder: "E.g. Pnma, Fd-3m, P6_3mc"})
registerFilter('results.material.symmetry.point_group', {...termQuantity, placeholder: "E.g. 6mm, m-3m, 6/mmm"})
registerFilter('results.material.symmetry.hall_symbol', {...termQuantity, placeholder: "E.g. F 4d 2 3 -1d"})
registerFilter('results.material.symmetry.prototype_aflow_id', {...termQuantity, placeholder: "E.g. A_cF8_227_a"})
registerFilter(
  'results.material.topology',
  nestedQuantity,
  [
    {name: 'label', ...termQuantity},
    {name: 'method', ...termQuantity},
    {name: 'description', ...termQuantity},
    {name: 'material_id', ...termQuantity},
    {name: 'material_name', ...termQuantity},
    {name: 'structural_type', ...termQuantity},
    {name: 'dimensionality', ...termQuantity},
    {name: 'building_block', ...termQuantity},
    {name: 'elements', ...termQuantity},
    {name: 'n_elements', ...numberHistogramQuantity},
    {name: 'chemical_formula_hill', ...termQuantity, placeholder: "E.g. H2O2, C2H5Br"},
    {name: 'chemical_formula_iupac', ...termQuantity, placeholder: "E.g. GaAs, SiC", label: 'Chemical formula IUPAC'},
    {name: 'chemical_formula_reduced', ...termQuantity, placeholder: "E.g. H2NaO, ClNa"},
    {name: 'chemical_formula_anonymous', ...termQuantity, placeholder: "E.g. A2B, A3B2C2"},
    {name: 'atomic_fraction', ...numberHistogramQuantity},
    {name: 'mass_fraction', ...numberHistogramQuantity},
    {name: 'n_atoms', ...numberHistogramQuantity},
    {name: 'system_relation.type', ...termQuantity},
    {name: 'cell.a', ...numberHistogramQuantity},
    {name: 'cell.b', ...numberHistogramQuantity},
    {name: 'cell.c', ...numberHistogramQuantity},
    {name: 'cell.alpha', ...numberHistogramQuantity},
    {name: 'cell.beta', ...numberHistogramQuantity},
    {name: 'cell.gamma', ...numberHistogramQuantity},
    {name: 'cell.atomic_density', ...numberHistogramQuantity},
    {name: 'cell.mass_density', ...numberHistogramQuantity},
    {name: 'cell.volume', ...numberHistogramQuantity},
    {name: 'symmetry.bravais_lattice', ...termQuantity},
    {name: 'symmetry.crystal_system', ...termQuantity},
    {name: 'symmetry.hall_symbol', ...termQuantity},
    {name: 'symmetry.point_group', ...termQuantity, placeholder: "E.g. 6mm, m-3m, 6/mmm"},
    {name: 'symmetry.space_group_symbol', ...termQuantity, placeholder: "E.g. Pnma, Fd-3m, P6_3mc"},
    {name: 'symmetry.strukturbericht_designation', ...termQuantity, placeholder: "E.g. Pnma, Fd-3m, P6_3mc"},
    {name: 'symmetry.prototype_name', ...termQuantity},
    {name: 'symmetry.prototype_label_aflow', ...termQuantity},
    {name: 'sbu_type', ...termQuantity},
    {name: 'largest_cavity_diameter', ...numberHistogramQuantity},
    {name: 'pore_limiting_diameter', ...numberHistogramQuantity},
    {name: 'largest_included_sphere_along_free_sphere_path', ...numberHistogramQuantity},
    {name: 'accessible_surface_area', ...numberHistogramQuantity},
    {name: 'accessible_volume', ...numberHistogramQuantity},
    {name: 'void_fraction', ...numberHistogramQuantity},
    {name: 'n_channels', ...numberHistogramQuantity},
    {name: 'sbu_coordination_number', ...numberHistogramQuantity}
  ]
)
registerFilter(
  'results.material.elemental_composition',
  nestedQuantity,
  [
    {name: 'element', ...termQuantity},
    {name: 'atomic_fraction', ...numberHistogramQuantity},
    {name: 'mass_fraction', ...numberHistogramQuantity}
  ]
)
registerFilter(
  'results.material.topology.elemental_composition',
  nestedQuantity,
  [
    {name: 'element', ...termQuantity},
    {name: 'mass', ...numberHistogramQuantity},
    {name: 'atomic_fraction', ...numberHistogramQuantity},
    {name: 'mass_fraction', ...numberHistogramQuantity}
  ]
)
registerFilter(
  'results.material.topology.active_orbitals',
  nestedQuantity,
  [
    {name: 'n_quantum_number', ...termQuantity},
    {name: 'l_quantum_number', ...termQuantity},
    {name: 'l_quantum_symbol', ...termQuantity},
    {name: 'ml_quantum_symbol', ...termQuantity},
    {name: 'ms_quantum_symbol', ...termQuantity},
    {name: 'j_quantum_number', ...termQuantity},
    {name: 'mj_quantum_number', ...termQuantity},
    {name: 'occupation', ...termQuantity},
    {name: 'n_electrons_excited', ...termQuantity},
    {name: 'degeneracy', ...termQuantity}
  ]
)
registerFilter('results.method.method_name', {...termQuantity, scale: 'log'})
registerFilter('results.method.workflow_name', {...termQuantity, scale: 'log'})
registerFilter('results.method.simulation.program_name', {...termQuantity, scale: 'log'})
registerFilter('results.method.simulation.program_version', termQuantity)
registerFilter('results.method.simulation.program_version_internal', termQuantity)
registerFilter('results.method.simulation.precision.native_tier', {...termQuantity, placeholder: "E.g. VASP - accurate", label: 'Code-specific tier'})
registerFilter('results.method.simulation.precision.k_line_density', {...numberHistogramQuantity, scale: 'log', label: 'k-line density'})
registerFilter('results.method.simulation.precision.basis_set', {...termQuantity, scale: 'log'})
registerFilter('results.method.simulation.precision.planewave_cutoff', {...numberHistogramQuantity, label: 'Plane-wave cutoff', scale: 'log'})
registerFilter('results.method.simulation.precision.apw_cutoff', {...numberHistogramQuantity, label: 'APW cutoff', scale: 'log'})
registerFilter('results.method.simulation.dft.core_electron_treatment', termQuantity)
registerFilter('results.method.simulation.dft.jacobs_ladder', {...termQuantity, scale: 'log', label: 'Jacob\'s ladder'})
registerFilter('results.method.simulation.dft.xc_functional_type', {
  ...termQuantity,
  scale: 'log',
  label: 'Jacob\'s ladder',
  options: {
    'LDA': {label: 'LDA'},
    'GGA': {label: 'GGA'},
    'meta-GGA': {label: 'Meta-GGA'},
    'hyper-GGA': {label: 'Hyper-GGA'},
    'hybrid': {label: 'Hybrid'}
  }
})
registerFilter('results.method.simulation.dft.xc_functional_names', {...termQuantityNonExclusive, scale: 'log', label: 'XC functional names'})
registerFilter('results.method.simulation.dft.exact_exchange_mixing_factor', {...numberHistogramQuantity, scale: 'log'})
registerFilter('results.method.simulation.dft.hubbard_kanamori_model.u_effective', {...numberHistogramQuantity, scale: 'log'})
registerFilter('results.method.simulation.dft.relativity_method', termQuantity)
registerFilter('results.method.simulation.tb.type', {...termQuantity, scale: 'log'})
registerFilter('results.method.simulation.tb.localization_type', {...termQuantity, scale: 'log'})
registerFilter('results.method.simulation.gw.type', {...termQuantity, label: 'GW type'})
registerFilter('results.method.simulation.gw.starting_point_type', {
  ...termQuantity,
  scale: 'log',
  options: {
    'LDA': {label: 'LDA'},
    'GGA': {label: 'GGA'},
    'meta-GGA': {label: 'Meta-GGA'},
    'hyper-GGA': {label: 'Hyper-GGA'},
    'hybrid': {label: 'Hybrid'},
    'HF': {label: 'HF'}
  }
})
registerFilter('results.method.simulation.gw.basis_set_type', {...termQuantity, scale: 'log'})
registerFilter('results.method.simulation.bse.type', termQuantity)
registerFilter('results.method.simulation.bse.solver', termQuantity)
registerFilter('results.method.simulation.bse.starting_point_type', {
  ...termQuantity,
  scale: 'log',
  options: {
    'LDA': {label: 'LDA'},
    'GGA': {label: 'GGA'},
    'meta-GGA': {label: 'Meta-GGA'},
    'hyper-GGA': {label: 'Hyper-GGA'},
    'hybrid': {label: 'Hybrid'},
    'HF': {label: 'HF'}
  }
})
registerFilter('results.method.simulation.bse.basis_set_type', {...termQuantity, scale: 'log'})
registerFilter('results.method.simulation.bse.gw_type', {...termQuantity, scale: 'log', label: `GW type`})
registerFilter('results.method.simulation.dmft.impurity_solver_type', {...termQuantity})
registerFilter('results.method.simulation.dmft.magnetic_state', {...termQuantity})
registerFilter('results.method.simulation.dmft.inverse_temperature', {...numberHistogramQuantity, scale: 'log'})
registerFilter('results.method.simulation.dmft.u', {...numberHistogramQuantity, scale: 'log'})
registerFilter('results.method.simulation.dmft.jh', {...numberHistogramQuantity, label: `JH`, scale: 'log'})
registerFilter('results.method.simulation.dmft.analytical_continuation', {...termQuantity})
registerFilter('results.eln.sections', termQuantity)
registerFilter('results.eln.tags', termQuantity)
registerFilter('results.eln.methods', termQuantity)
registerFilter('results.eln.instruments', termQuantity)
registerFilter('results.eln.lab_ids', {...termQuantity, label: 'Lab IDs'})
registerFilter('results.eln.names', noAggQuantity)
registerFilter('results.eln.descriptions', noAggQuantity)
registerFilter('external_db', {...termQuantity, label: 'External database', scale: 'log'})
registerFilter('authors.name', {...termQuantityNonExclusive, label: 'Author name'})
registerFilter('upload_create_time', {...numberHistogramQuantity, scale: 'log'})
registerFilter('entry_create_time', {...numberHistogramQuantity, scale: 'log'})
registerFilter('datasets.dataset_name', {...termQuantityLarge, label: 'Dataset name'})
registerFilter('datasets.doi', {...termQuantity, label: 'Dataset DOI'})
registerFilter('datasets.dataset_id', termQuantity)
registerFilter('domain', termQuantity)
registerFilter('entry_id', termQuantity)
registerFilter('entry_name', termQuantity)
registerFilter('mainfile', termQuantity)
registerFilter('upload_id', termQuantity)
registerFilter('upload_name', termQuantity)
registerFilter('published', termQuantity)
registerFilter('main_author.user_id', termQuantity)
registerFilter('quantities', {...noAggQuantity, label: 'Metainfo definition', queryMode: 'all'})
registerFilter('sections', {...noAggQuantity, label: 'Metainfo sections', queryMode: 'all'})
registerFilter('files', {...noAggQuantity, queryMode: 'all'})
registerFilter('section_defs.definition_qualified_name', {...noAggQuantity, label: 'Section defs qualified name', queryMode: 'all'})
registerFilter('entry_references.target_entry_id', {...noAggQuantity, label: 'Entry references target entry id', queryMode: 'all'})
registerFilter('entry_type', {...noAggQuantity, label: 'Entry type', queryMode: 'all'})
registerFilter('entry_name.prefix', {...noAggQuantity, label: 'Entry name', queryMode: 'all'})
registerFilter('results.material.material_id', termQuantity)
registerFilter('optimade_filter', {multiple: true, queryMode: 'all'})
registerFilter('processed', {label: 'Processed', queryMode: 'all'})
registerFilter('text_search_contents', {multiple: true, queryMode: 'all'})
registerFilter('custom_quantities', {
  serializerExact: value => {
    const jsonStr = JSON.stringify(value)
    const result = encodeURIComponent(jsonStr)
    return result
  },
  serializerPretty: value => {
    if (!value) {
      return null
    }
    return `${value.and.length} condition${value.and.length > 0 ? 's' : ''}`
  },
  deserializer: value => {
    const jsonStr = decodeURIComponent(value)
    const result = JSON.parse(jsonStr)
    return result
  },
  multiple: false,
  value: {
    set: () => {
      // TODO: We ignore the query here, it is later added to the final API query.
      // We had to do this hack, because there is not way to add a logical query
      // behind a prefix.
    }
  }
})
registerFilter(
  'results.properties.spectroscopic.spectra.provenance.eels',
  {...nestedQuantity, label: 'Electron energy loss spectrum (EELS)'},
  [
    {name: 'detector_type', ...termQuantity},
    {name: 'resolution', ...numberHistogramQuantity},
    {name: 'min_energy', ...numberHistogramQuantity},
    {name: 'max_energy', ...numberHistogramQuantity}
  ]
)
registerFilter(
  'results.properties.electronic.band_structure_electronic',
  {...nestedQuantity, label: 'Band structure'},
  [
    {name: 'spin_polarized', label: 'Spin-polarized', ...termQuantityBool}
  ]
)
registerFilter(
  'results.properties.electronic.dos_electronic',
  {...nestedQuantity, label: 'Density of states'},
  [
    {name: 'spin_polarized', label: 'Spin-polarized', ...termQuantityBool}
  ]
)
registerFilter(
  'results.properties.electronic.band_structure_electronic.band_gap',
  nestedQuantity,
  [
    {name: 'type', ...termQuantity},
    {name: 'value', ...numberHistogramQuantity, scale: 'log'}
  ]
)
registerFilter(
  'results.properties.electronic.band_gap',
  nestedQuantity,
  [
    {name: 'type', ...termQuantity},
    {name: 'value', ...numberHistogramQuantity, scale: 'log'}
  ]
)
registerFilter('results.properties.electronic.band_gap.provenance.label', termQuantity)
registerFilter(
  'results.properties.optoelectronic.solar_cell',
  nestedQuantity,
  [
    {name: 'efficiency', ...numberHistogramQuantity, scale: 'log'},
    {name: 'fill_factor', ...numberHistogramQuantity, scale: 'log'},
    {name: 'open_circuit_voltage', ...numberHistogramQuantity, scale: 'log'},
    {name: 'short_circuit_current_density', ...numberHistogramQuantity, scale: 'log'},
    {name: 'illumination_intensity', ...numberHistogramQuantity, scale: 'log'},
    {name: 'device_area', ...numberHistogramQuantity, scale: 'log'},
    {name: 'device_architecture', ...termQuantity},
    {name: 'absorber_fabrication', ...termQuantity},
    {name: 'device_stack', ...termQuantityAllNonExclusive},
    {name: 'absorber', ...termQuantityAllNonExclusive},
    {name: 'electron_transport_layer', ...termQuantityAllNonExclusive},
    {name: 'hole_transport_layer', ...termQuantityAllNonExclusive},
    {name: 'substrate', ...termQuantityAllNonExclusive},
    {name: 'back_contact', ...termQuantityAllNonExclusive}
  ]
)
registerFilter(
  'results.properties.catalytic.catalyst',
  nestedQuantity,
  [
    {name: 'characterization_methods', ...termQuantity},
    {name: 'surface_area', ...numberHistogramQuantity, scale: 'log'},
    {name: 'catalyst_name', ...termQuantity},
    {name: 'catalyst_type', ...termQuantity},
    {name: 'preparation_method', ...termQuantity}
  ]
)
registerFilter(
  'results.properties.catalytic.reaction',
  nestedQuantity,
  [
    {name: 'name', ...termQuantity},
    {name: 'type', ...termQuantity}
  ]
)
registerFilter(
  'results.properties.catalytic.reaction.reaction_conditions',
  nestedQuantity,
  [
    {name: 'temperature', ...numberHistogramQuantity, scale: 'linear'},
    {name: 'pressure', ...numberHistogramQuantity, scale: 'linear'},
    {name: 'weight_hourly_space_velocity', ...numberHistogramQuantity, scale: 'linear'},
    {name: 'gas_hourly_space_velocity', ...numberHistogramQuantity, scale: 'linear'},
    {name: 'flow_rate', ...numberHistogramQuantity, scale: 'linear'},
    {name: 'time_on_stream', ...numberHistogramQuantity, scale: 'linear'}
  ]
)
registerFilter(
  'results.properties.catalytic.reaction.products',
  nestedQuantity,
  [
    {name: 'name', ...termQuantityAllNonExclusive},
    {name: 'gas_concentration_out', ...numberHistogramQuantity, scale: 'linear'},
    {name: 'selectivity', ...numberHistogramQuantity, scale: 'linear'},
    {name: 'space_time_yield', ...numberHistogramQuantity, scale: 'linear'}
  ]
)
registerFilter(
  'results.properties.catalytic.reaction.reactants',
  nestedQuantity,
  [
    {name: 'name', ...termQuantityAllNonExclusive},
    {name: 'gas_concentration_in', ...numberHistogramQuantity, scale: 'linear'},
    {name: 'gas_concentration_out', ...numberHistogramQuantity, scale: 'linear'},
    {name: 'conversion', ...numberHistogramQuantity, scale: 'linear'}
  ]
)
registerFilter(
  'results.properties.catalytic.reaction.rates',
  nestedQuantity,
  [
    {name: 'name', ...termQuantityAllNonExclusive},
    {name: 'reaction_rate', ...numberHistogramQuantity, scale: 'linear'},
    {name: 'specific_mass_rate', ...numberHistogramQuantity, scale: 'linear'},
    {name: 'specific_surface_area_rate', ...numberHistogramQuantity, scale: 'linear'},
    {name: 'rate', ...numberHistogramQuantity, scale: 'linear'},
    {name: 'turnover_frequency', ...numberHistogramQuantity, scale: 'linear'}
  ]
)
registerFilter(
  'results.properties.catalytic.reaction.reaction_mechanism',
  nestedQuantity,
  [
    {name: 'initial_states', ...termQuantity},
    {name: 'final_states', ...termQuantity},
    {name: 'reaction_enthalpy', ...numberHistogramQuantity, scale: 'linear'},
    {name: 'activation_energy', ...numberHistogramQuantity, scale: 'linear'}
  ]
)
registerFilter(
  'results.properties.mechanical.bulk_modulus',
  nestedQuantity,
  [
    {name: 'type', ...termQuantity},
    {name: 'value', ...numberHistogramQuantity}
  ]
)
registerFilter(
  'results.properties.mechanical.shear_modulus',
  nestedQuantity,
  [
    {name: 'type', ...termQuantity},
    {name: 'value', ...numberHistogramQuantity}
  ]
)
registerFilter(
  'results.properties.available_properties',
  termQuantityAll
)
registerFilter(
  'results.properties.mechanical.energy_volume_curve',
  nestedQuantity,
  [
    {name: 'type', ...termQuantity}
  ]
)
registerFilter(
  'results.properties.geometry_optimization',
  nestedQuantity,
  [
    {name: 'final_energy_difference', ...numberHistogramQuantity, scale: 'log'},
    {name: 'final_displacement_maximum', ...numberHistogramQuantity, scale: 'log'},
    {name: 'final_force_maximum', ...numberHistogramQuantity, scale: 'log'}
  ]
)
registerFilter(
  'results.properties.thermodynamic.trajectory',
  nestedQuantity,
  [
    {name: 'available_properties', ...termQuantityAll},
    {name: 'provenance.molecular_dynamics.ensemble_type', ...termQuantity},
    {name: 'provenance.molecular_dynamics.time_step', ...numberHistogramQuantity}
  ]
)

// Visibility: controls the 'owner'-parameter in the API query, not part of the
// query itself.
registerFilter(
  'visibility',
  {
    ...noQueryQuantity,
    default: 'visible',
    description: 'The visibility of the entry.'
  }
)

// Combine: controls whether materials search combines data from several
// entries.
registerFilter(
  'combine',
  {
    ...noQueryQuantity,
    default: true,
    description: 'If selected, your filters may be matched from several entries that contain the same material. When unchecked, the material has to have a single entry that matches all your filters.'
  }
)

// Exclusive: controls the way elements search is done.
registerFilter(
  'exclusive',
  {
    ...noQueryQuantity,
    dtype: DType.Boolean,
    default: false,
    description: "Search for entries with compositions that only (exclusively) contain the selected atoms. The default is to return all entries that have at least (inclusively) the selected atoms."
  }
)

// In exclusive element query the elements names are sorted and concatenated
// into a single string.
registerFilter(
  'results.material.elements',
  {
    widget: {
      search_quantity: 'results.material.elements',
      type: 'periodic_table',
      scale: 'log',
      layout: {
        sm: {w: 12, h: 8, minW: 12, minH: 8},
        md: {w: 12, h: 8, minW: 12, minH: 8},
        lg: {w: 12, h: 8, minW: 12, minH: 8},
        xl: {w: 12, h: 8, minW: 12, minH: 8},
        xxl: {w: 12, h: 8, minW: 12, minH: 8}
      }
    },
    aggs: {terms: {size: elementData.elements.length}},
    value: {
      set: (newQuery, oldQuery, value) => {
        if (oldQuery.exclusive) {
          if (value.size !== 0) {
            newQuery['results.material.elements_exclusive'] = setToArray(value).sort().join(' ')
          }
        } else {
          newQuery['results.material.elements'] = value
        }
      }
    },
    multiple: true,
    exclusive: false,
    queryMode: 'all'
  }
)

// Electronic properties: subset of results.properties.available_properties
registerFilterOptions(
  'electronic_properties',
  'results.properties.available_properties',
  'Electronic properties',
  'The electronic properties that are present in an entry.',
  {
    'electronic.band_structure_electronic.band_gap': {label: 'Band gap'},
    band_structure_electronic: {label: 'Band structure'},
    dos_electronic: {label: 'Density of states'},
    greens_functions_electronic: {label: 'Green\u0027s functions'},
    eels: {label: 'Electron energy loss spectrum'}
  }
)

// Vibrational properties: subset of results.properties.available_properties
registerFilterOptions(
  'vibrational_properties',
  'results.properties.available_properties',
  'Vibrational properties',
  'The vibrational properties that are present in an entry.',
  {
    dos_phonon: {label: 'Phonon density of states'},
    band_structure_phonon: {label: 'Phonon band structure'},
    energy_free_helmholtz: {label: 'Helmholtz free energy'},
    heat_capacity_constant_volume: {label: 'Heat capacity constant volume'}
  }
)

// Mechanical properties: subset of results.properties.available_properties
registerFilterOptions(
  'mechanical_properties',
  'results.properties.available_properties',
  'Mechanical properties',
  'The mechanical properties that are present in an entry.',
  {
    bulk_modulus: {label: 'Bulk modulus'},
    shear_modulus: {label: 'Shear modulus'},
    energy_volume_curve: {label: 'Energy-volume curve'}
  }
)

// Spectroscopic properties: subset of results.properties.available_properties
registerFilterOptions(
  'spectroscopic_properties',
  'results.properties.available_properties',
  'Spectroscopic properties',
  'The spectroscopic properties that are present in an entry.',
  {
    eels: {label: 'Electron energy loss spectrum'}
  }
)

// Thermodynamical properties: subset of results.properties.available_properties
registerFilterOptions(
  'thermodynamic_properties',
  'results.properties.available_properties',
  'Thermodynamic properties',
  'The thermodynamic properties that are present.',
  {
    trajectory: {label: 'Trajectory'}
  }
)

/**
 * Creates static suggestion for all metainfo quantities that have an enum
 * value. Also provides suggestions for quantity names.
 */
export const quantityNameSearch = 'quantity name'
export function getStaticSuggestions(quantities, filterData) {
  const suggestions = {}
  const filters = quantities
    ? [...quantities]
    : Object.keys(filterData)

  // Add suggestions from metainfo
  for (const quantity of filters) {
    const data = searchQuantities[quantity]
    const isEnum = data?.type?.type_kind?.toLowerCase() === 'enum'
    if (isEnum) {
      const options = data.type.type_data
      const maxLength = Math.max(...options.map(option => option.length))
      const minLength = maxLength <= 2 ? 1 : 2
      suggestions[quantity] = getSuggestions(
        options,
        minLength
      )
    }
  }

  // Add suggestions for quantity names
  if (quantities.has(quantityNameSearch)) {
    suggestions[quantityNameSearch] = getSuggestions(
      filters.filter(value => value !== quantityNameSearch && !filterData[value].section),
      2
    )
  }
  return suggestions
}

/**
 * HOC that is used to preload search quantities from all required schemas. This
 * simplifies the rendering logic by first loading all schemas before rendering
 * any components that rely on them.
 */
export const withSearchQuantities = (WrappedComponent) => {
  const WithFilters = ({initialSearchQuantities, ...rest}) => {
    // Here we load the python schemas, and determine which YAML/Nexus schemas
    // to load later.
    const [yamlOptions, nexusOptions, initialFilterData] = useMemo(() => {
      const options = getOptions(initialSearchQuantities)
      const yamlOptions = options.filter((name) => name.includes(`#${yamlSchemaPrefix}`))
      const nexusOptions = options.filter((name) => name.startsWith('nexus.'))

      // Perform glob filtering on default filters. Only exclude affects the
      // default filters.
      const defaultFilters = {}
      for (const [key, value] of Object.entries(defaultFilterData)) {
        if (glob(key, [key], initialSearchQuantities?.exclude)) {
          defaultFilters[key] = value
        }
      }

      // Load the python filters from plugins
      const pythonFilterData = {}
      for (const [name, def] of Object.entries(searchQuantities)) {
        if (def.dynamic && glob(name, initialSearchQuantities?.include, initialSearchQuantities?.exclude)) {
          const {path, schema} = parseQuantityName(name)
          const params = {
            name: path,
            quantity: name,
            schema,
            aggregatable: def.aggregatable
          }
          pythonFilterData[name] = new Filter(def, params)
        }
      }

      return [
        yamlOptions,
        nexusOptions,
        {...defaultFilters, ...pythonFilterData}
      ]
    }, [initialSearchQuantities])

    const metainfo = useGlobalMetainfo()
    const [loadingYaml, setLoadingYaml] = useState(yamlOptions.length)
    const [loadingNexus, setLoadingNexus] = useState(nexusOptions.length)
    const [filters, setFilters] = useState(initialFilterData)
    const { raiseError } = useErrors()

    // Nexus metainfo is loaded here once metainfo is ready
    useEffect(() => {
      if (!nexusOptions.length || !metainfo) return
      const pkg = metainfo._packageDefs['nexus']
      const sections = pkg.section_definitions
      const nexusFilters = {}
      for (const section of sections) {
        const sectionPath = `nexus.${section.name}`

        // The NeXus section is skipped (it contains duplicate information)
        if (sectionPath === 'nexus.NeXus') continue

        // Only applications definitions are loaded
        if (section?.more?.nx_category !== 'application') continue

        // Sections from which no quantities are included are skipped
        if (!glob(sectionPath, initialSearchQuantities?.include, initialSearchQuantities?.exclude) && !initialSearchQuantities?.include.some(x => x.includes(sectionPath))) {
          continue
        }

        // Add all included quantities recursively
        for (const [def, path, repeats] of getQuantities(section)) {
          const filterPath = `${sectionPath}.${path}`
          const included = glob(filterPath, initialSearchQuantities?.include, initialSearchQuantities?.exclude)
          if (!included) continue
          const dtype = dtypeMap[getDatatype(def)]
          // TODO: For some Nexus quantities, the data types cannot be fetched.
          if (!dtype) {
            continue
          }
          const params = {
            name: filterPath,
            quantity: filterPath,
            aggregatable: new Set([DType.String, DType.Enum, DType.Boolean]).has(getDatatype(def)),
            repeats: repeats
          }
          nexusFilters[filterPath] = new Filter(def, params)
        }
      }
      setFilters((old) => ({...old, ...nexusFilters}))
      setLoadingNexus(false)
    }, [metainfo, nexusOptions, initialSearchQuantities])

    // YAML schemas are loaded here asynchronously
    useEffect(() => {
      if (!yamlOptions.length || !metainfo) return
      async function fetchSchemas(options) {
        const yamlFilters = {}
        for (const schemaPath of options) {
          let schemaDefinition
          try {
            schemaDefinition = await metainfo.resolveDefinition(schemaPath)
          } catch (e) {
            raiseError(`
              Unable to load the schema ${schemaPath} that is used in this app.
              If the schema has not been published, please make sure that you
              are logged in and have the correct access rights.
            `)
            throw e
          }
          for (const [def, path] of getQuantities(schemaDefinition)) {
            const dtype = dtypeMap[getDatatype(def)]
            if (!dtype) {
              throw Error(`Unable to load the data type for ${path}.`)
            }
            const filterPath = `data.${path}${schemaSeparator}${schemaPath}`
            const included = glob(filterPath, initialSearchQuantities?.include, initialSearchQuantities?.exclude)
            if (!included) continue
            const apiPath = `${filterPath}${dtypeSeparator}${dtype}`
            const {path: quantity, schema} = parseQuantityName(filterPath)
            yamlFilters[filterPath] = new Filter(def, {name: path, schema, quantity, requestQuantity: apiPath})
          }
        }
        setFilters((old) => ({...old, ...yamlFilters}))
        setLoadingYaml(false)
      }

      // Get a list of all distinct YAML schemas. These are loaded one by one
      // and the required filters are registered from each.
      const yamlSchemas = new Set(yamlOptions.map(x => x.split('#').pop()))
      fetchSchemas(yamlSchemas)
    }, [yamlOptions, metainfo, raiseError, initialSearchQuantities])

    return (loadingYaml || loadingNexus)
      ? <Box margin={1}>
          <Typography>Loading the required schemas...</Typography>
        </Box>
      : <WrappedComponent {...rest} initialSearchQuantities={filters}/>
  }

  WithFilters.displayName = `withFilter(${WrappedComponent.displayName || WrappedComponent.name})`
  WithFilters.propTypes = {
    initialSearchQuantities: PropTypes.object // Determines which filters are available
  }

  return WithFilters
}
