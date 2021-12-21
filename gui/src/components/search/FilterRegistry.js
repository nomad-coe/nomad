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
import { isNil, isString } from 'lodash'
import { setToArray, getDatatype, getSerializer, getDeserializer } from '../../utils'
import searchQuantities from '../../searchQuantities'
import { getDimension } from '../../units'
import InputList from './input/InputList'
import InputPeriodicTable from './input/InputPeriodicTable'
import elementData from '../../elementData'

// Containers for filter information
export const filterGroups = [] // Mapping from a group name -> set of filter names
export const filterAbbreviations = [] // Mapping of filter full name -> abbreviation
export const filterFullnames = [] // Mapping of filter abbreviation -> full name
export const filterDataGlobal = {} // Stores data for each registered filter

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
export const labelVibrational = 'Vibrational'
export const labelMechanical = 'Mechanical'
export const labelSpectroscopy = 'Spectroscopy'
export const labelAuthor = 'Author / Origin'
export const labelAccess = 'Access'
export const labelDataset = 'Dataset'
export const labelIDs = 'IDs'
export const labelArchive = 'Archive'

/**
 * This function is used to register a new filter within the SearchContext.
 * Filters are entities that can be searched through the filter panel and the
 * search bar, and can be encoded in the URL. Notice that a filter in this
 * context does not have to correspond to a quantity in the metainfo.
 *
 * Only registered filters may be searched for. The registration must happen
 * before any components use the filters. This is because:
 *  - The initial aggregation results must be fetched before any components
 *  using the filter values are rendered.
 *  - Several components need to know the list of available filters (e.g. the
 *  search bar and  the search panel). If filters are only registered during
 *  component initialization, it may already be too late to update other
 *  components.
 *
 * @param {string} name Name of the filter.
 * @param {string} group The group into which the filter belongs to. Groups
 * are used to e.g. in showing FilterSummaries about a group of filters.
 * @param {quantity} options Data object containing options for the filter. Can
 * include the following data:
  *  - agg: Object containing a custom setter/getter for the aggregation value.
  *      As a shortcut you can provide an ES aggregation type as a string,
  *  - value: Object containfig a custom setter/getter for the filter value.
  *  - multiple: Whether the user can simultaneously provide multiple values for
  *      this filter.
  *  - exclusive: Whether this filter is exclusive: only one value may be
  *      associated with an entry.
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
  *  - queryMode: The default query mode (e.g. 'any', 'all) when multiple values
  *    can specified for this filter. Defaults to 'any'.
  *  - guiOnly: Whether this filter is only shown in the GUI and does not get
  *    serialized into the API call.
  *  - default: A default value which is implicitly enforced in the API call.
  *    This value will not be serialized in the search bar.
  *  - resources: A list of resources for which this filter is enabled.
  *  - aggSize: Aggregation size.
  */
function registerFilter(name, group, quantity, subQuantities) {
  function save(name, group, quantity) {
    if (group) {
      filterGroups[group]
        ? filterGroups[group].add(name)
        : filterGroups[group] = new Set([name])
    }

    const data = filterDataGlobal[name] || {}
    const agg = quantity.agg
    if (agg) {
      let aggSet, aggGet
      if (isString(agg)) {
        aggSet = {[name]: {type: agg}}
        aggGet = (aggs) => (aggs[name][agg].data)
      } else {
        aggSet = agg.set
        aggGet = agg.get
      }
      data.aggSet = aggSet
      data.aggGet = aggGet
    }
    if (quantity.value) {
      data.valueSet = quantity.value.set
    }
    data.multiple = quantity.multiple === undefined ? true : quantity.multiple
    data.exclusive = quantity.exclusive === undefined ? true : quantity.exclusive
    data.stats = quantity.stats
    data.options = quantity.options
    data.unit = quantity.unit || searchQuantities[name]?.unit
    data.dtype = quantity.dtype || getDatatype(name)
    data.serializerExact = getSerializer(data.dtype, false)
    data.serializerPretty = getSerializer(data.dtype, true)
    data.dimension = getDimension(data.unit)
    data.deserializer = getDeserializer(data.dtype, data.dimension)
    data.label = quantity.label
    data.section = !isNil(searchQuantities[name]?.nested)
    data.description = quantity.description
    data.scale = quantity.scale || 1
    data.aggSize = quantity.aggSize
    data.aggSizeOverride = quantity.aggSizeOverride
    if (data.queryMode && !data.multiple) {
      throw Error('Only filters that accept multiple values may have a query mode.')
    }
    data.queryMode = quantity.queryMode || 'any'
    data.guiOnly = quantity.guiOnly
    if (quantity.default && !data.guiOnly) {
      throw Error('Only filters that do not correspond to a metainfo value may have default values set.')
    }
    data.default = quantity.default
    data.resources = new Set(quantity.resources || ['entries', 'materials'])
    filterDataGlobal[name] = data
  }
  save(name, group, quantity)

  // Register section subquantities
  if (subQuantities) {
    for (let quantity of subQuantities) {
      let subname = `${name}.${quantity.name}`
      save(subname, group, quantity)
    }
  }
}

// Configuration for the docked statistics
const listStatConfig = {
  component: InputList,
  layout: {
    width: 'small',
    ratio: 3 / 4
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
const termQuantity = {agg: 'terms', stats: listStatConfig, aggSize: 5}
const termQuantityBool = {
  agg: 'terms',
  stats: listStatConfig,
  aggSize: 2,
  options: {
    false: {label: 'false'},
    true: {label: 'true'}
  }
}
const termQuantityNonExclusive = {agg: 'terms', stats: listStatConfig, exclusive: false, aggSize: 5}
const noAggQuantity = {stats: listStatConfig}
const nestedQuantity = {}
const noQueryQuantity = {guiOnly: true, multiple: false}
const rangeQuantity = {agg: 'min_max', multiple: false}

// Filters that directly correspond to a metainfo value
registerFilter('results.material.structural_type', labelMaterial, {...termQuantity, scale: 1 / 4})
registerFilter('results.material.functional_type', labelMaterial, termQuantityNonExclusive)
registerFilter('results.material.compound_type', labelMaterial, termQuantityNonExclusive)
registerFilter('results.material.material_name', labelMaterial, termQuantity)
registerFilter('results.material.chemical_formula_hill', labelElements, termQuantity)
registerFilter('results.material.chemical_formula_anonymous', labelElements, termQuantity)
registerFilter('results.material.n_elements', labelElements, {...rangeQuantity, label: 'Number of Elements'})
registerFilter('results.material.symmetry.bravais_lattice', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.crystal_system', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.structure_name', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.strukturbericht_designation', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.space_group_symbol', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.point_group', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.hall_symbol', labelSymmetry, termQuantity)
registerFilter('results.material.symmetry.prototype_aflow_id', labelSymmetry, termQuantity)
registerFilter('results.method.method_name', labelMethod, {...termQuantity, scale: 1 / 4})
registerFilter('results.method.simulation.program_name', labelSimulation, {...termQuantity, scale: 1 / 4})
registerFilter('results.method.simulation.program_version', labelSimulation, termQuantity)
registerFilter('results.method.simulation.dft.basis_set_type', labelDFT, {...termQuantity, scale: 1 / 4})
registerFilter('results.method.simulation.dft.core_electron_treatment', labelDFT, termQuantity)
registerFilter('results.method.simulation.dft.xc_functional_type', labelDFT, {...termQuantity, scale: 1 / 2, label: 'XC Functional Type'})
registerFilter('results.method.simulation.dft.relativity_method', labelDFT, termQuantity)
registerFilter('results.method.simulation.gw.type', labelGW, {...termQuantity, label: 'GW Type'})
registerFilter('external_db', labelAuthor, {...termQuantity, label: 'External Database'})
registerFilter('authors.name', labelAuthor, {...termQuantity, label: 'Author Name'})
registerFilter('upload_create_time', labelAuthor, rangeQuantity)
registerFilter('datasets.dataset_name', labelDataset, {...termQuantity, label: 'Dataset Name'})
registerFilter('datasets.doi', labelDataset, {...noAggQuantity, label: 'Dataset DOI'})
registerFilter('entry_id', labelIDs, noAggQuantity)
registerFilter('upload_id', labelIDs, noAggQuantity)
registerFilter('quantities', labelArchive, {...noAggQuantity, label: 'Metainfo definition', queryMode: 'all'})
registerFilter('results.material.material_id', labelIDs, noAggQuantity)
registerFilter('datasets.dataset_id', labelIDs, noAggQuantity)
registerFilter(
  'results.properties.spectroscopy.eels',
  labelSpectroscopy,
  {...nestedQuantity, label: 'Electron Energy Loss Spectrum (EELS)'},
  [
    {name: 'detector_type', ...termQuantity},
    {name: 'resolution', ...rangeQuantity},
    // {name: 'min_energy', ...rangeQuantity},
    // {name: 'max_energy', ...rangeQuantity}
    {
      name: 'energy_window',
      stats: listStatConfig,
      agg: {
        set: {
          'results.properties.spectroscopy.eels.min_energy': {
            type: 'min_max',
            exclude: (updated) => updated?.has('results.properties.spectroscopy.eels.energy_window')
          },
          'results.properties.spectroscopy.eels.max_energy': {
            type: 'min_max',
            exclude: (updated) => updated?.has('results.properties.spectroscopy.eels.energy_window')
          }
        },
        get: (aggs) => {
          const min = aggs['results.properties.spectroscopy.eels.min_energy']
          const max = aggs['results.properties.spectroscopy.eels.max_energy']
          return [min.min_max.data[0], max.min_max.data[1]]
        }
      },
      value: {
        set: (newQuery, oldQuery, value) => {
          newQuery['results.properties.spectroscopy.eels'] = {
            min_energy: {
              gte: value.gte
            },
            max_energy: {
              lte: value.lte
            }
          }
        }
      },
      multiple: false,
      exlusive: true,
      dtype: 'number',
      unit: 'joule',
      label: 'Energy Window',
      description: 'Defines bounds for the minimum and maximum energies in the spectrum.'
    }
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
    {name: 'value', ...rangeQuantity}
  ]
)
registerFilter(
  'results.properties.mechanical.bulk_modulus',
  labelMechanical,
  nestedQuantity,
  [
    {name: 'type', ...termQuantity},
    {name: 'value', ...rangeQuantity}
  ]
)
registerFilter(
  'results.properties.mechanical.shear_modulus',
  labelMechanical,
  nestedQuantity,
  [
    {name: 'type', ...termQuantity},
    {name: 'value', ...rangeQuantity}
  ]
)
registerFilter(
  'results.properties.mechanical.energy_volume_curve',
  labelMechanical,
  nestedQuantity,
  [
    {name: 'type', ...termQuantity}
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
    agg: 'terms',
    aggSize: elementData.elements.length,
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
const electronicOptions = {
  'electronic.band_structure_electronic.band_gap': {label: 'Band gap'},
  band_structure_electronic: {label: 'Band structure'},
  dos_electronic: {label: 'Density of states'}
}
const electronicProps = new Set(Object.keys(electronicOptions))
registerFilter(
  'electronic_properties',
  labelElectronic,
  {
    stats: listStatConfig,
    agg: {
      set: {'results.properties.available_properties': {type: 'terms'}},
      get: (aggs) => (aggs['results.properties.available_properties'].terms.data
        .filter((value) => electronicProps.has(value.value)))
    },
    aggSizeOverride: 200,
    value: {
      set: (newQuery, oldQuery, value) => {
        const data = newQuery['results.properties.available_properties'] || new Set()
        value.forEach((item) => { data.add(item) })
        newQuery['results.properties.available_properties:all'] = data
      }
    },
    multiple: true,
    exclusive: false,
    options: electronicOptions,
    label: 'Electronic properties',
    description: 'The electronic properties that are present in an entry.'
  }
)
// Vibrational properties: subset of results.properties.available_properties
export const vibrationalOptions = {
  dos_phonon: {label: 'Phonon density of states'},
  band_structure_phonon: {label: 'Phonon band structure'},
  energy_free_helmholtz: {label: 'Helmholtz free energy'},
  heat_capacity_constant_volume: {label: 'Heat capacity constant volume'}
}
const vibrationalProps = new Set(Object.keys(vibrationalOptions))
registerFilter(
  'vibrational_properties',
  labelVibrational,
  {
    stats: listStatConfig,
    agg: {
      set: {'results.properties.available_properties': {type: 'terms'}},
      get: (aggs) => (aggs['results.properties.available_properties'].terms.data
        .filter((value) => vibrationalProps.has(value.value)))
    },
    aggSizeOverride: 200,
    value: {
      set: (newQuery, oldQuery, value) => {
        const data = newQuery['results.properties.available_properties'] || new Set()
        value.forEach((item) => { data.add(item) })
        newQuery['results.properties.available_properties:all'] = data
      }
    },
    multiple: true,
    exclusive: false,
    options: vibrationalOptions,
    label: 'Vibrational properties',
    description: 'The vibrational properties that are present in an entry.'
  }
)
// Mechanical properties: subset of results.properties.available_properties
export const mechanicalOptions = {
  bulk_modulus: {label: 'Bulk modulus'},
  shear_modulus: {label: 'Shear modulus'},
  energy_volume_curve: {label: 'Energy-volume curve'}
}
const mechanicalProps = new Set(Object.keys(mechanicalOptions))
registerFilter(
  'mechanical_properties',
  labelMechanical,
  {
    stats: listStatConfig,
    agg: {
      set: {'results.properties.available_properties': {type: 'terms'}},
      get: (aggs) => (aggs['results.properties.available_properties'].terms.data
        .filter((value) => mechanicalProps.has(value.value)))
    },
    aggSizeOverride: 200,
    value: {
      set: (newQuery, oldQuery, value) => {
        const data = newQuery['results.properties.available_properties'] || new Set()
        value.forEach((item) => { data.add(item) })
        newQuery['results.properties.available_properties:all'] = data
      }
    },
    multiple: true,
    exclusive: false,
    options: mechanicalOptions,
    label: 'Mechanical properties',
    description: 'The mechanical properties that are present in an entry.'
  }
)
// Spectroscopic properties: subset of results.properties.available_properties
export const spectroscopicOptions = {
  eels: {label: 'Electron energy loss spectrum'}
}
const spectroscopicProps = new Set(Object.keys(spectroscopicOptions))
registerFilter(
  'spectroscopic_properties',
  labelSpectroscopy,
  {
    stats: listStatConfig,
    agg: {
      set: {'results.properties.available_properties': {type: 'terms'}},
      get: (aggs) => (aggs['results.properties.available_properties'].terms.data
        .filter((value) => spectroscopicProps.has(value.value)))
    },
    aggSizeOverride: 200,
    value: {
      set: (newQuery, oldQuery, value) => {
        const data = newQuery['results.properties.available_properties'] || new Set()
        value.forEach((item) => { data.add(item) })
        newQuery['results.properties.available_properties:all'] = data
      }
    },
    multiple: true,
    exclusive: false,
    options: spectroscopicOptions,
    label: 'Spectroscopic properties',
    description: 'The spectroscopic properties that are present in an entry.'
  }
)
// The filter abbreviation mapping has to be done only after all filters have
// been registered.
const abbreviations = {}
const nameAbbreviationPairs = [...Object.keys(filterDataGlobal)].map(
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
