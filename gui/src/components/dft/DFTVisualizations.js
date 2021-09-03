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
import { resolveRef } from '../archive/metainfo'
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
      'dft.workflow.type'
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
      <Grid item xs={12}>
        <QuantityHistogram quantity="dft.labels_springer_classification" title="Functional classification" initialScale={1} multiple/>
      </Grid>
      <Grid item xs={12}>
        <QuantityHistogram quantity="dft.workflow.type" title="Workflows" valueLabels={workflowTypeLabels} initialScale={0.25} />
      </Grid>
    </Grid>
  )
}

DFTPropertyVisualizations.propTypes = {
  info: PropTypes.object
}
