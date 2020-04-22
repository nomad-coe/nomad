import React, { useContext, useEffect } from 'react'
import PropTypes from 'prop-types'
import { Grid } from '@material-ui/core'
import { Quantity } from '../search/QuantityHistogram'
import { searchContext } from '../search/SearchContext'

export function DFTMethodVisualizations(props) {
  const {info} = props
  const {response: {statistics, metric}, setStatistics} = useContext(searchContext)
  useEffect(() => {
    setStatistics(['dft.code_name', 'dft.basis_set', 'dft.xc_functional'])
  }, [])

  if (statistics.code_name && info) {
    // filter based on known codes, since elastic search might return 0 aggregations on
    // obsolete code names
    const filteredCodeNames = {}
    const defaultValue = {
      code_runs: 0
    }
    defaultValue[metric] = 0
    info.codes.forEach(key => {
      filteredCodeNames[key] = statistics.code_name[key] || defaultValue
    })
    statistics.code_name = filteredCodeNames
  }

  return (
    <Grid container spacing={2}>
      <Grid item xs={8}>
        <Quantity quantity="dft.code_name" title="Code" scale={0.25} metric={metric} sort columns={2} />
      </Grid>
      <Grid item xs={4}>
        <Quantity quantity="dft.basis_set" title="Basis set" scale={0.25} metric={metric} sort />
        <Quantity quantity="dft.xc_functional" title="XC functionals" scale={0.5} metric={metric} sort />
      </Grid>
    </Grid>
  )
}

DFTMethodVisualizations.propTypes = {
  info: PropTypes.object
}

export function DFTSystemVisualizations(props) {
  const {info} = props
  const {response: {statistics, metric}, setStatisticsToRefresh, setStatistics} = useContext(searchContext)
  useEffect(() => {
    setStatisticsToRefresh('dft.labels_springer_compound_class')
    setStatistics(['dft.labels_springer_compound_class', 'dft.system', 'dft.crystal_system', 'dft.compound_type'])
  }, [])

  if (statistics.code_name && info) {
    // filter based on known codes, since elastic search might return 0 aggregations on
    // obsolete code names
    const filteredCodeNames = {}
    const defaultValue = {
      code_runs: 0
    }
    defaultValue[metric] = 0
    info.codes.forEach(key => {
      filteredCodeNames[key] = statistics.code_name[key] || defaultValue
    })
    statistics.code_name = filteredCodeNames
  }

  return (
    <Grid container spacing={2}>
      <Grid item xs={4}>
        <Quantity quantity="dft.compound_type" title="Compound type" scale={1} metric={metric} sort />
      </Grid>
      <Grid item xs={4}>
        <Quantity quantity="dft.system" title="System type" scale={0.25} metric={metric} sort />
        <Quantity quantity="dft.crystal_system" title="Crystal system" scale={1} metric={metric} sort />
      </Grid>
      <Grid item xs={4}>
        <Quantity quantity="dft.labels_springer_compound_class" title="Springer compound" scale={1} metric={metric} />
      </Grid>
    </Grid>
  )
}

DFTSystemVisualizations.propTypes = {
  info: PropTypes.object
}

const searchable_quantities_categories = {
  energy_quantities: [
    'energy_total',
    'energy_total_T0',
    'energy_free',
    'energy_electrostatic',
    'energy_X',
    'energy_XC',
    'energy_sum_eigenvalues'
  ],
  electronic_quantities: [
    'dos_values',
    'eigenvalues_values',
    'volumetric_data_values',
    'electronic_kinetic_energy',
    'total_charge',
    'atomic_multipole_values'
  ],
  forces_quantities: [
    'atom_forces_free',
    'atom_forces_raw',
    'atom_forces_T0',
    'atom_forces',
    'stress_tensor'
  ],
  vibrational_quantities: [
    'thermodynamical_property_heat_capacity_C_v',
    'vibrational_free_energy_at_constant_volume',
    'band_energies'
  ],
  magnetic_quantities: [
    'spin_S2'
  ],
  optical_quantities: [
    'excitation_energies',
    'oscillator_strengths',
    'transition_dipole_moments'
  ]
}

export function DFTPropertyVisualizations(props) {
  const {info} = props
  const {response: {statistics, metric}, setStatisticsToRefresh, setStatistics} = useContext(searchContext)
  useEffect(() => {
    setStatisticsToRefresh('dft.labels_springer_classification')
    setStatistics([
      'dft.searchable_quantities',
      'dft.labels_springer_classification'
    ])
  }, [])

  if (statistics.code_name && info) {
    // filter based on known codes, since elastic search might return 0 aggregations on
    // obsolete code names
    const filteredCodeNames = {}
    const defaultValue = {
      code_runs: 0
    }
    defaultValue[metric] = 0
    info.codes.forEach(key => {
      filteredCodeNames[key] = statistics.code_name[key] || defaultValue
    })
    statistics.code_name = filteredCodeNames
  }

  const data = (category) => {
    const results = {}
    const data = statistics['dft.searchable_quantities']
    if (!data) {
      return null
    }

    searchable_quantities_categories[category].forEach(value => {
      if (data[value]) {
        results[value] = data[value]
      }
    })
    return results
  }

  return (
    <Grid container spacing={2}>
      <Grid item xs={4}>
        <Quantity quantity="dft.searchable_quantities" data={data('energy_quantities')} title="Energy" scale={1} metric={metric} sort tooltips />
        <Quantity quantity="dft.searchable_quantities" data={data('electronic_quantities')} title="Electronic" scale={1} metric={metric} sort tooltips />
      </Grid>
      <Grid item xs={4}>
        <Quantity quantity="dft.searchable_quantities" data={data('forces_quantities')} title="Forces" scale={1} metric={metric} sort tooltips />
        <Quantity quantity="dft.searchable_quantities" data={data('vibrational_quantities')} title="Vibrational" scale={1} metric={metric} sort tooltips />
        <Quantity quantity="dft.searchable_quantities" data={data('optical_quantities')} title="Optical" scale={1} metric={metric} sort tooltips />
      </Grid>
      <Grid item xs={4}>
        <Quantity quantity="dft.labels_springer_classification" title="Springer classification" scale={1} metric={metric} tooltips />
        <Quantity quantity="dft.searchable_quantities" data={data('magnetic_quantities')} title="Magnetic" scale={1} metric={metric} sort tooltips />
      </Grid>
    </Grid>
  )
}

DFTPropertyVisualizations.propTypes = {
  info: PropTypes.object
}
