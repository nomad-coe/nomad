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
import { isNil, isArray } from 'lodash'
import { setToArray, getDatatype, getSerializer, getDeserializer, formatLabel, DType } from '../../utils'
import searchQuantities from '../../searchQuantities'
import { Unit } from '../../units'
import InputList from './input/InputList'
import InputRange from './input/InputRange'
import InputSection from './input/InputSection'
import InputPeriodicTable from './input/InputPeriodicTable'
import elementData from '../../elementData'

// Containers for filter information
export const filterGroups = [] // Mapping from a group name -> set of filter names
export const filterAbbreviations = [] // Mapping of filter full name -> abbreviation
export const filterFullnames = [] // Mapping of filter abbreviation -> full name
export const filterData = {} // Stores data for each registered filter

// Labels for the filter menus
export const labelMaterial = 'Material'
export const labelElements = 'Elements / Formula'
export const labelSymmetry = 'Symmetry'
export const labelMethod = 'Method'
export const labelSimulation = 'Simulation'
export const labelDFT = 'DFT'
export const labelGW = 'GW'
export const labelExperiment = 'Experiment'
export const labelEELS = 'EELS'
export const labelProperties = 'Properties'
export const labelElectronic = 'Electronic'
export const labelOptoelectronic = 'Optoelectronic'
export const labelVibrational = 'Vibrational'
export const labelMechanical = 'Mechanical'
export const labelSpectroscopy = 'Spectroscopy'
export const labelThermodynamic = 'Thermodynamic'
export const labelGeometryOptimization = 'Geometry optimization'
export const labelELN = 'Electronic Lab Notebook'
export const labelAuthor = 'Author / Origin'
export const labelAccess = 'Access'
export const labelDataset = 'Dataset'
export const labelIDs = 'IDs'
export const labelArchive = 'Processed data quantities'

/**
 * Used to gather a list of fixed filter options from the metainfo.
 * @param {string} quantity Metainfo name
 * @returns Dictionary containing the available options and their labels.
 */
function getEnumOptions(quantity) {
  const metainfoOptions = searchQuantities?.[quantity]?.type?.type_data
  if (isArray(metainfoOptions) && metainfoOptions.length > 0) {
    const opt = {}
    for (const name of metainfoOptions) {
      opt[name] = {label: name}
    }

    // We do not display the option for 'not processed': it is more of a
    // debug value
    delete opt['not processed']
    return opt
  }
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
 * @param {string} group The group into which the filter belongs to. Groups
 * are used to e.g. in showing FilterSummaries about a group of filters.
 * @param {obj} config Data object containing options for the filter. Can
 * include the following data:
  *  - aggs: Object containing default values for specific aggregation types
  *      that may be requested for this filter.
  *      Also completely customized setter/getter methods are supported. E.g.
  *        aggs = {
  *          terms: {size: 5},
  *          histogram: {buckets: 20},
  *          min_max: {set: (config) => ({}), get: (agg) => ({})}
  *        }
  *  - value: Object containing a custom setter/getter for the filter value.
  *  - multiple: Whether the user can simultaneously provide multiple values for
  *      this filter.
  *  - queryMode: The default query mode (e.g. 'any', 'all) when multiple values
  *      can specified for this filter. Defaults to 'any'.
  *  - exclusive: Whether this filter is exclusive: only one value is
  *      associated with a single entry.
  *  - stats: Object that determines how this filter is visualized if it is
  *      docked above the search results.
  *  - options: Object containing explicit options that this filter supports.
  *  - unit: The unit for this filter. If no value is given and the name
  *      corresponds to a metainfo name the data type is read directly from the
  *      metainfo.
  *  - dtype: The data type for this filter. If no value is given and the
  *      name corresponds to a metainfo name the data type is read directly from
  *      the metainfo.
  *  - label: Name of the filter shown in the GUI. If no value is given and the
  *      name corresponds to a metainfo name the description is read directly
  *      from the metainfo.
  *  - description: Description of the filter shown e.g. in the tooltips. If no
  *      value is given and the name corresponnds to a metainfo, the metainfo
  *      description is used.
  *  - guiOnly: Whether this filter is only shown in the GUI and does not get
  *      serialized into the API call.
  *  - default: A default value which is implicitly enforced in the API call.
  *      This value will not be serialized in the search bar.
  *  - resources: A list of resources (entries, materials) for which this filter
  *      is enabled. Default is that a filter is available both for entries and
  *      materials.
  */
function saveFilter(name, group, config, parent) {
  if (group) {
    filterGroups[group]
      ? filterGroups[group].add(name)
      : filterGroups[group] = new Set([name])
  }

  const data = filterData[name] || {}
  data.options = config.options || getEnumOptions(name)
  data.aggs = config.aggs
  data.value = config.value
  data.multiple = config.multiple === undefined ? true : config.multiple
  data.exclusive = config.exclusive === undefined ? true : config.exclusive
  data.unit = config.unit || searchQuantities[name]?.unit
  data.dtype = config.dtype || getDatatype(name)
  data.serializerExact = getSerializer(data.dtype, false)
  data.serializerPretty = getSerializer(data.dtype, true)
  data.dimension = data.unit && new Unit(data.unit).dimension()
  data.deserializer = getDeserializer(data.dtype, data.dimension)
  data.label = config.label || formatLabel(searchQuantities[name]?.name || name)
  data.nested = searchQuantities[name]?.nested
  data.section = !isNil(data.nested)
  data.repeats = searchQuantities[name]?.repeats
  data.stats = config.stats || getStatsComponent(parent, data.dtype)
  data.description = config.description || searchQuantities[name]?.description
  data.scale = config.scale || 'linear'
  if (data.queryMode && !data.multiple) {
    throw Error('Only filters that accept multiple values may have a query mode.')
  }
  data.queryMode = config.queryMode || 'any'
  data.guiOnly = config.guiOnly
  if (config.default && !data.guiOnly) {
    throw Error('Only filters that do not correspond to a metainfo value may have default values set.')
  }
  data.default = config.default
  data.resources = new Set(config.resources || ['entries', 'materials'])
  filterData[name] = data
}

/**
 * Used to register a filter value for an individual metainfo quantity or
 * section.
 */
function registerFilter(name, group, config, subQuantities) {
  saveFilter(name, group, config)
  if (subQuantities) {
    for (const subConfig of subQuantities) {
      const subname = `${name}.${subConfig.name}`
      saveFilter(subname, group, subConfig, name)
    }
  }
}

/**
 * Used to register a filter that is based on a subset of quantity values.
 */
function registerFilterOptions(name, group, target, label, description, options) {
  const keys = Object.keys(options)
  registerFilter(
    name,
    group,
    {
      aggs: {
        terms: {
          set: ((target) => (config) => ({quantity: target, include: keys, ...config}))(target),
          get: ((target) => (agg) => agg)(target)
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
      exclusive: false,
      queryMode: 'all',
      options: options,
      label: label,
      description: description
    }
  )
}

/**
 * Tries to automatically create a suitable statistics component for the given
 * quantity.
 * @param {string} quantity
 * @param {bool} nested
 * @param {DType} dtype
 * @returns A statistics setup.
 */
const getStatsComponent = (parent, dtype) => {
  const section = filterData?.[parent]?.section
  const wrap = (InputComp) => {
    return section
      ? (props) => {
        return <InputSection
          section={parent}
          visible
          disableHeader
        >
          <InputComp {...props}/>
        </InputSection>
      }
      : InputComp
  }
  if (dtype === DType.Float || dtype === DType.Int || dtype === DType.Timestamp) {
    return {
      component: wrap(InputRange),
      layout: {
        width: 'medium',
        ratio: 3 / 1
      }
    }
  } else {
    return {
      component: wrap(InputList),
      layout: {
        width: 'small',
        ratio: 3 / 4
      }
    }
  }
}

const ptStatConfig = {
  component: InputPeriodicTable,
  layout: {
    width: 'large',
    ratio: 3 / 2
  }
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
const noAggQuantity = {}
const nestedQuantity = {}
const noQueryQuantity = {guiOnly: true, multiple: false}
const numberHistogramQuantity = {multiple: false, exclusive: false}

// Filters that directly correspond to a metainfo value
registerFilter('results.material.structural_type', labelMaterial, {...termQuantity, scale: '1/4'})
registerFilter('results.material.functional_type', labelMaterial, termQuantityNonExclusive)
registerFilter('results.material.compound_type', labelMaterial, termQuantityNonExclusive)
registerFilter('results.material.material_name', labelMaterial, termQuantity)
registerFilter('results.material.chemical_formula_hill', labelElements, termQuantity)
registerFilter('results.material.chemical_formula_anonymous', labelElements, termQuantity)
registerFilter('results.material.n_elements', labelElements, {...numberHistogramQuantity, label: 'Number of Elements'})
registerFilter('results.material.symmetry.bravais_lattice', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.crystal_system', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.structure_name', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.strukturbericht_designation', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.space_group_symbol', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.point_group', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.hall_symbol', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.prototype_aflow_id', labelSymmetry, termQuantity)
registerFilter('results.method.method_name', labelMethod, {...termQuantity, scale: '1/4'})
registerFilter('results.method.workflow_name', labelMethod, {...termQuantity, scale: '1/4'})
registerFilter('results.method.simulation.program_name', labelSimulation, {...termQuantity, scale: '1/4'})
registerFilter('results.method.simulation.program_version', labelSimulation, termQuantity)
registerFilter('results.method.simulation.dft.basis_set_type', labelDFT, {...termQuantity, scale: '1/4'})
registerFilter('results.method.simulation.dft.core_electron_treatment', labelDFT, termQuantity)
registerFilter('results.method.simulation.dft.xc_functional_type', labelDFT, {...termQuantity, scale: '1/2', label: 'XC Functional Type'})
registerFilter('results.method.simulation.dft.xc_functional_names', labelDFT, {...termQuantityNonExclusive, scale: '1/2', label: 'XC Functional Names'})
registerFilter('results.method.simulation.dft.relativity_method', labelDFT, termQuantity)
registerFilter('results.method.simulation.gw.type', labelGW, {...termQuantity, label: 'GW Type'})
registerFilter('results.eln.sections', labelELN, termQuantity)
registerFilter('results.eln.tags', labelELN, termQuantity)
registerFilter('results.eln.methods', labelELN, termQuantity)
registerFilter('results.eln.instruments', labelELN, termQuantity)
registerFilter('results.eln.lab_ids', labelELN, termQuantity)
registerFilter('results.eln.names', labelELN, termQuantity)
registerFilter('results.eln.descriptions', labelELN, termQuantity)
registerFilter('external_db', labelAuthor, {...termQuantity, label: 'External Database', scale: '1/4'})
registerFilter('authors.name', labelAuthor, {...termQuantityNonExclusive, label: 'Author Name'})
registerFilter('upload_create_time', labelAuthor, {...numberHistogramQuantity, scale: '1/2'})
registerFilter('datasets.dataset_name', labelDataset, {...termQuantityLarge, label: 'Dataset Name'})
registerFilter('datasets.doi', labelDataset, {...termQuantity, label: 'Dataset DOI'})
registerFilter('entry_id', labelIDs, termQuantity)
registerFilter('upload_id', labelIDs, termQuantity)
registerFilter('quantities', labelArchive, {...noAggQuantity, label: 'Metainfo definition', queryMode: 'all'})
registerFilter('results.material.material_id', labelIDs, termQuantity)
registerFilter('datasets.dataset_id', labelIDs, termQuantity)
registerFilter(
  'results.properties.spectroscopy.eels',
  labelSpectroscopy,
  {...nestedQuantity, label: 'Electron Energy Loss Spectrum (EELS)'},
  [
    {name: 'detector_type', ...termQuantity},
    {name: 'resolution', ...numberHistogramQuantity},
    {name: 'min_energy', ...numberHistogramQuantity},
    {name: 'max_energy', ...numberHistogramQuantity}
  ]
)
registerFilter(
  'results.properties.electronic.band_structure_electronic',
  labelElectronic,
  {...nestedQuantity, label: 'Band Structure'},
  [
    {name: 'spin_polarized', label: 'Spin-polarized', ...termQuantityBool}
  ]
)
registerFilter(
  'results.properties.electronic.dos_electronic',
  labelElectronic,
  {...nestedQuantity, label: 'Density of States (DOS)'},
  [
    {name: 'spin_polarized', label: 'Spin-polarized', ...termQuantityBool}
  ]
)
registerFilter(
  'results.properties.electronic.band_structure_electronic.band_gap',
  labelElectronic,
  nestedQuantity,
  [
    {name: 'type', ...termQuantity},
    {name: 'value', ...numberHistogramQuantity, scale: '1/4'}
  ]
)
registerFilter(
  'results.properties.optoelectronic.band_gap_optical',
  labelOptoelectronic,
  nestedQuantity,
  [
    {name: 'type', ...termQuantity},
    {name: 'value', ...numberHistogramQuantity, scale: '1/4'}
  ]
)
registerFilter(
  'results.properties.optoelectronic.solar_cell',
  labelOptoelectronic,
  nestedQuantity,
  [
    {name: 'efficiency', ...numberHistogramQuantity, scale: '1/4'},
    {name: 'fill_factor', ...numberHistogramQuantity, scale: '1/4'},
    {name: 'open_circuit_voltage', ...numberHistogramQuantity, scale: '1/4'},
    {name: 'short_circuit_current_density', ...numberHistogramQuantity, scale: '1/4'},
    {name: 'illumination_intensity', ...numberHistogramQuantity, scale: '1/4'},
    {name: 'device_area', ...numberHistogramQuantity, scale: '1/4'},
    {name: 'device_architecture', ...termQuantity},
    {name: 'device_stack', ...termQuantity},
    {name: 'absorber', ...termQuantity},
    {name: 'absorber_fabrication', ...termQuantity},
    {name: 'electron_transport_layer', ...termQuantity},
    {name: 'hole_transport_layer', ...termQuantity},
    {name: 'substrate', ...termQuantity},
    {name: 'back_contact', ...termQuantity}
  ]
)
registerFilter(
  'results.properties.mechanical.bulk_modulus',
  labelMechanical,
  nestedQuantity,
  [
    {name: 'type', ...termQuantity},
    {name: 'value', ...numberHistogramQuantity}
  ]
)
registerFilter(
  'results.properties.mechanical.shear_modulus',
  labelMechanical,
  nestedQuantity,
  [
    {name: 'type', ...termQuantity},
    {name: 'value', ...numberHistogramQuantity}
  ]
)
registerFilter(
  'results.properties.available_properties',
  labelProperties,
  termQuantityAll
)
registerFilter(
  'results.properties.mechanical.energy_volume_curve',
  labelMechanical,
  nestedQuantity,
  [
    {name: 'type', ...termQuantity}
  ]
)
registerFilter(
  'results.properties.geometry_optimization',
  labelGeometryOptimization,
  nestedQuantity,
  [
    {name: 'final_energy_difference', ...numberHistogramQuantity, scale: '1/8'},
    {name: 'final_displacement_maximum', ...numberHistogramQuantity, scale: '1/8'},
    {name: 'final_force_maximum', ...numberHistogramQuantity, scale: '1/8'}
  ]
)
registerFilter(
  'results.properties.thermodynamic.trajectory',
  labelThermodynamic,
  nestedQuantity,
  [
    {name: 'available_properties', ...termQuantityAll},
    {name: 'methodology.molecular_dynamics.ensemble_type', ...termQuantity},
    {name: 'methodology.molecular_dynamics.time_step', ...numberHistogramQuantity}
  ]
)

// Visibility: controls the 'owner'-parameter in the API query, not part of the
// query itself.
registerFilter('visibility', labelAccess, {...noQueryQuantity, default: 'visible'})

// Combine: controls whether materials search combines data from several
// entries.
registerFilter('combine', undefined, {
  ...noQueryQuantity,
  default: true,
  resources: ['materials']
})

// Exclusive: controls the way elements search is done.
registerFilter('exclusive', undefined, {...noQueryQuantity, default: false})

// In exclusive element query the elements names are sorted and concatenated
// into a single string.
registerFilter(
  'results.material.elements',
  labelElements,
  {
    stats: ptStatConfig,
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
  labelElectronic,
  'results.properties.available_properties',
  'Electronic Properties',
  'The electronic properties that are present in an entry.',
  {
    'electronic.band_structure_electronic.band_gap': {label: 'Band gap'},
    band_structure_electronic: {label: 'Band structure'},
    dos_electronic: {label: 'Density of states'}
  }
)

// Optoelectronic properties: subset of results.properties.available_properties
registerFilterOptions(
  'optoelectronic_properties',
  labelOptoelectronic,
  'results.properties.available_properties',
  'Optoelectronic properties',
  'The optoelectronic properties that are present in an entry.',
  {
    'band_gap_optical': {label: 'Band gap optical'},
    'solar_cell': {label: 'Solar cell'}
  }
)

// Vibrational properties: subset of results.properties.available_properties
registerFilterOptions(
  'vibrational_properties',
  labelVibrational,
  'results.properties.available_properties',
  'Vibrational Properties',
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
  labelMechanical,
  'results.properties.available_properties',
  'Mechanical Properties',
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
  labelSpectroscopy,
  'results.properties.available_properties',
  'Spectroscopic Properties',
  'The spectroscopic properties that are present in an entry.',
  {
    eels: {label: 'Electron energy loss spectrum'}
  }
)

// Thermodynamical properties: subset of results.properties.available_properties
registerFilterOptions(
  'thermodynamic_properties',
  labelThermodynamic,
  'results.properties.available_properties',
  'Thermodynamic Properties',
  'The thermodynamic properties that are present.',
  {
    trajectory: {label: 'Trajectory'}
  }
)

// The filter abbreviation mapping has to be done only after all filters have
// been registered.
const abbreviations = {}
const nameAbbreviationPairs = [...Object.keys(filterData)].map(
  fullname => [fullname, fullname.split('.').pop()])
for (const [fullname, abbreviation] of nameAbbreviationPairs) {
  const old = abbreviations[abbreviation]
  if (old === undefined) {
    abbreviations[abbreviation] = 1
  } else {
    abbreviations[abbreviation] += 1
  }
  filterAbbreviations[fullname] = fullname
  filterFullnames[fullname] = fullname
}
for (const [fullname, abbreviation] of nameAbbreviationPairs) {
  if (abbreviations[abbreviation] === 1) {
    filterAbbreviations[fullname] = abbreviation
    filterFullnames[abbreviation] = fullname
  }
}

// Material and entry queries target slightly different fields. Here we prebuild
// the mapping.
export const materialNames = {} // Mapping of field name from entry -> material
export const entryNames = {} // Mapping of field name from material -> entry
for (const name of Object.keys(searchQuantities)) {
  const prefix = 'results.material.'
  let materialName
  if (name.startsWith(prefix)) {
    materialName = name.substring(prefix.length)
  } else {
    materialName = `entries.${name}`
  }
  materialNames[name] = materialName
  entryNames[materialName] = name
}

/**
 * Function for creating static suggestions. Mimics the suggestion logic used by
 * the suggestions API endpoint.
 *
 * @param {str} category Category for the suggestions
 * @param {array} values Array of available values
 * @param {number} minLength Minimum input length before suggestions are considered.
 * @param {func} text Function that maps the value into the suggested text input
 *
 * @return {object} Object containing a list of options and a function for
 *   filtering them based on the input.
 */
function getSuggestions(
  category, values, minLength = 2, text = (value) => value) {
  const options = values
    .map(value => {
      const optionCleaned = value.trim().replace(/_/g, ' ').toLowerCase()
      const matches = [...optionCleaned.matchAll(/[ .]/g)]
      let tokens = [optionCleaned]
      tokens = tokens.concat(matches.map(match => optionCleaned.slice(match.index + 1)))
      return {
        value: value,
        category: category,
        text: text && text(value),
        tokens: tokens
      }
    })
  const filter = (input) => {
    // Minimum input length
    if (input.length < minLength) {
      return []
    }
    // Gather all matches
    const inputCleaned = input.trim().replace(/_/g, ' ').toLowerCase()
    let suggestions = options.filter(option => option.tokens.some(token => token.startsWith(inputCleaned)))

    // Sort matches based on value length (the more the input covers from the
    // value, the better the match)
    suggestions = suggestions.sort((a, b) => a.value.length - b.value.length)
    return suggestions
  }

  return {options, filter}
}

/**
 * Creates static suggestion for all metainfo quantities that have an enum
 * value. Also provides suggestions for quantity names.
 */
function getStaticSuggestions() {
  const suggestions = {}
  const filters = Object.keys(filterData)

  // Add suggestions from metainfo
  for (const quantity of filters) {
    const data = searchQuantities[quantity]
    const isEnum = data?.type?.type_kind === 'Enum'
    if (isEnum) {
      const options = data.type.type_data
      const maxLength = Math.max(...options.map(option => option.length))
      const minLength = maxLength <= 2 ? 1 : 2
      suggestions[quantity] = getSuggestions(
        quantity,
        options,
        minLength,
        (value) => `${quantity}=${value}`
      )
    }
  }

  // Add suggestions for quantity names
  suggestions['quantity name'] = getSuggestions(
    'quantity name',
    filters.filter(value => !filterData[value].section),
    2
  )
  return suggestions
}

export const staticSuggestions = getStaticSuggestions()
