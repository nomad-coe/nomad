import React, { useContext, useEffect } from 'react'
import PropTypes from 'prop-types'
import { Grid } from '@material-ui/core'
import QuantityHistogram from '../search/QuantityHistogram'
import TooltipQuantityHistogram from '../search/TooltipQuantityHistogram'
import { searchContext } from '../search/SearchContext'

export function DFTMethodVisualizations(props) {
  const {info} = props
  const {response: {statistics, metric}, setStatistics} = useContext(searchContext)
  useEffect(() => {
    setStatistics(['dft.code_name', 'dft.basis_set', 'dft.xc_functional'])
    // eslint-disable-next-line
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
        <QuantityHistogram quantity="dft.code_name" title="Code" initialScale={0.25} columns={2} />
      </Grid>
      <Grid item xs={4}>
        <QuantityHistogram quantity="dft.basis_set" title="Basis set" initialScale={0.25} />
        <QuantityHistogram quantity="dft.xc_functional" title="XC functionals" initialScale={0.5} />
      </Grid>
    </Grid>
  )
}

DFTMethodVisualizations.propTypes = {
  info: PropTypes.object
}

export function DFTSystemVisualizations(props) {
  const {info} = props
  const {response: {statistics, metric}, setStatistics} = useContext(searchContext)
  useEffect(() => {
    setStatistics(['dft.labels_springer_compound_class', 'dft.system', 'dft.crystal_system', 'dft.compound_type'])
    // eslint-disable-next-line
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
        <QuantityHistogram quantity="dft.system" title="System type" initialScale={0.25} />
        <QuantityHistogram quantity="dft.crystal_system" title="Crystal system" />
      </Grid>
      <Grid item xs={4}>
        <QuantityHistogram quantity="dft.compound_type" title="Compound type" initialScale={0.25} />
      </Grid>
      <Grid item xs={4}>
        <QuantityHistogram quantity="dft.labels_springer_compound_class" title="Compound classification" />
      </Grid>
    </Grid>
  )
}

DFTSystemVisualizations.propTypes = {
  info: PropTypes.object
}

const electronic_quantities = [
  'electronic_band_structure',
  'electronic_dos',
  'eigenvalues_values',
]
const mechanical_quantities = [
  'stress_tensor'
]
const thermal_quantities = [
  'thermodynamical_property_heat_capacity_C_v',
  'vibrational_free_energy_at_constant_volume',
  'phonon_band_structure',
  'phonon_dos'
]
const magnetic_quantities = [
  'spin_S2'
]
const optical_quantities = [
  'oscillator_strengths',
  'transition_dipole_moments'
]

const labels = {
  'energy_total': 'Total energy',
  'energy_total_T0': 'Total energy (0K)',
  'energy_free': 'Free energy',
  'energy_electrostatic': 'Electrostatic',
  'energy_X': 'Exchange',
  'energy_XC': 'Exchange-correlation',
  'energy_sum_eigenvalues': 'Band energy',
  'dos_values': 'DOS',
  'eigenvalues_values': 'Eigenvalues',
  'volumetric_data_values': 'Volumetric data',
  'electronic_kinetic_energy': 'Kinetic energy',
  'total_charge': 'Charge',
  'atom_forces_free': 'Free atomic forces',
  'atom_forces_raw': 'Raw atomic forces',
  'atom_forces_T0': 'Atomic forces (0K)',
  'atom_forces': 'Atomic forces',
  'stress_tensor': 'Stress tensor',
  'electronic_band_structure': 'Electronic band structure',
  'electronic_dos': 'Electronic density of states',
  'phonon_band_structure': 'Phonon dispersion',
  'phonon_dos': 'Phonon density of states',
  'thermodynamical_property_heat_capacity_C_v': 'Heat capacity',
  'vibrational_free_energy_at_constant_volume': 'Helmholtz free energy',
  'band_energies': 'Band energies',
  'spin_S2': 'Angular spin momentum squared',
  'excitation_energies': 'Excitation energies',
  'oscillator_strengths': 'Oscillator strengths',
  'transition_dipole_moments': 'Transition dipole moments',
  'atomic_multipole_values': 'Atomic multipole values'}

const workflowTypeLabels = {
  'geometry_optimization': 'Geometry optimization',
  'phonon': 'Phonons',
  'elastic': 'Elastic constants'
}

export function DFTPropertyVisualizations(props) {
  const {info} = props
  const {response: {statistics, metric}, setStatistics} = useContext(searchContext)
  useEffect(() => {
    setStatistics([
      'dft.searchable_quantities',
      'dft.labels_springer_classification',
      'dft.workflow.workflow_type'
    ])
    // eslint-disable-next-line
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
      <Grid item xs={7}>
        <QuantityHistogram quantity="dft.searchable_quantities" values={electronic_quantities} valueLabels={labels} title="Electronic" initialScale={0.5} tooltips multiple/>
        <QuantityHistogram quantity="dft.searchable_quantities" values={mechanical_quantities} valueLabels={labels} title="Mechanical" initialScale={0.5} multiple/>
        <QuantityHistogram quantity="dft.searchable_quantities" values={thermal_quantities} valueLabels={labels} title="Thermal" initialScale={0.5} multiple/>
        <QuantityHistogram quantity="dft.searchable_quantities" values={optical_quantities} valueLabels={labels} title="Optical" initialScale={1} multiple/>
        <QuantityHistogram quantity="dft.searchable_quantities" values={magnetic_quantities} valueLabels={labels} title="Magnetic" initialScale={1} multiple/>
      </Grid>
      <Grid item xs={5}>
        <QuantityHistogram quantity="dft.labels_springer_classification" title="Functional classification" initialScale={1} multiple/>
      </Grid>
      <Grid item xs={12}>
        <QuantityHistogram quantity="dft.workflow.workflow_type" title="Workflows" valueLabels={workflowTypeLabels} initialScale={0.25} />
      </Grid>
    </Grid>
  )
}

DFTPropertyVisualizations.propTypes = {
  info: PropTypes.object
}
