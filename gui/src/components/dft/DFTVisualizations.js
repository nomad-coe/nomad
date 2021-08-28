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
import React, { useContext, useEffect } from 'react'
import PropTypes from 'prop-types'
import { Grid } from '@material-ui/core'
import { makeStyles } from '@material-ui/core/styles'
import QuantityHistogram from '../search/QuantityHistogram'
import { searchContext } from '../search/SearchContext'
import { path, defsByName, resolveRef } from '../archive/metainfo'
import { nomadTheme } from '../../config'
import Markdown from '../Markdown'

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
  'eigenvalues_values'
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
  'eigenvalues_values': 'eigenvalues',
  'stress_tensor': 'stress tensor',
  'electronic_band_structure': 'electronic band structure',
  'electronic_dos': 'electronic density of states',
  'phonon_band_structure': 'phonon dispersion',
  'phonon_dos': 'phonon density of states',
  'thermodynamical_property_heat_capacity_C_v': 'heat capacity',
  'vibrational_free_energy_at_constant_volume': 'helmholtz free energy',
  'spin_S2': 'angular spin momentum squared',
  'oscillator_strengths': 'oscillator strengths',
  'transition_dipole_moments': 'transition dipole moments'
}
const metainfoPaths = {
  'eigenvalues_values': path('eigenvalues_values'),
  'stress_tensor': path('stress_tensor'),
  'electronic_band_structure': 'EntryArchive/section_run/section_single_configuration_calculation/section_k_band',
  'electronic_dos': 'EntryArchive/section_run/section_single_configuration_calculation/section_dos',
  'phonon_band_structure': 'EntryArchive/section_run/section_single_configuration_calculation/section_k_band',
  'phonon_dos': 'EntryArchive/section_run/section_single_configuration_calculation/section_dos',
  'thermodynamical_property_heat_capacity_C_v': path('thermodynamical_property_heat_capacity_C_v'),
  'vibrational_free_energy_at_constant_volume': path('vibrational_free_energy_at_constant_volume'),
  'spin_S2': path('spin_S2'),
  'oscillator_strengths': path('oscillator_strengths'),
  'transition_dipole_moments': path('transition_dipole_moments')
}

const useMetainfoTooltipStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    flexDirection: 'column',
    padding: 2
  },
  tooltipMarkdown: {
    fontSize: nomadTheme.overrides.MuiTooltip.tooltip.fontSize,
    color: 'white',
    '& a': {
      color: theme.palette.primary.light
    }
  }
}))

function MetaInfoTooltip({def, path}) {
  const classes = useMetainfoTooltipStyles()
  let description = def.description
  if (!description && def.sub_section) {
    description = resolveRef(def.sub_section)?.description
  }
  return <div className={classes.root} >
    <Markdown
      classes={{root: classes.tooltipMarkdown}}
    >{`${description?.slice(0, description.indexOf('.') || undefined)}. Click [here](/metainfo/${path}) for full the definition.`}</Markdown>
  </div>
}

MetaInfoTooltip.propTypes = {
  def: PropTypes.object,
  path: PropTypes.string
}

const tooltips = {}
for (const label in labels) {
  const path = metainfoPaths[label]
  console.assert(path, `The metainfo for ${label} does not exist.`)
  if (path) {
    const realName = path.split('/').slice(-1)[0]
    const metainfoDef = defsByName[realName]?.[0]
    console.assert(metainfoDef, `No metainfo definition called ${realName} does exist.`)
    if (metainfoDef) {
      tooltips[label] = <MetaInfoTooltip def={metainfoDef} path={path}></MetaInfoTooltip>
    }
  }
}

const workflowTypeLabels = {
  'geometry_optimization': 'geometry optimization',
  'phonon': 'phonons',
  'elastic': 'elastic constants',
  'molecular_dynamics': 'molecular dynamics'
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
        <QuantityHistogram quantity="dft.searchable_quantities" values={electronic_quantities} valueLabels={labels} title="Electronic" initialScale={0.5} tooltips={tooltips} multiple/>
        <QuantityHistogram quantity="dft.searchable_quantities" values={mechanical_quantities} valueLabels={labels} title="Mechanical" initialScale={0.5} tooltips={tooltips} multiple/>
        <QuantityHistogram quantity="dft.searchable_quantities" values={thermal_quantities} valueLabels={labels} title="Thermal" initialScale={0.5} tooltips={tooltips} multiple/>
        <QuantityHistogram quantity="dft.searchable_quantities" values={optical_quantities} valueLabels={labels} title="Optical" initialScale={1} tooltips={tooltips} multiple/>
        <QuantityHistogram quantity="dft.searchable_quantities" values={magnetic_quantities} valueLabels={labels} title="Magnetic" initialScale={1} tooltips={tooltips} multiple/>
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
