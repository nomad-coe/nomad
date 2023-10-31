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
import React, {useState, useEffect} from 'react'
import PropTypes from 'prop-types'
import { useTheme, makeStyles} from '@material-ui/core/styles'
import { Box } from '@material-ui/core'
import Plot from '../plotting/Plot'
import {
    PropertyProvenanceItem,
    PropertyProvenanceList
} from '../entry/properties/PropertyCard'
import { withErrorHandler } from '../ErrorHandler'
import { mergeObjects, alphabet, resolveInternalRef } from '../../utils'

export const gfError = 'Could not load Greens functions.'
const types = {
  greens_function_iw_im: {
    xLabel: 'i\u03c9<sub>n</sub>',
    yLabel: 'Im G(i\u03c9<sub>n</sub>)',
    axis: 'matsubara_freq'
  },
  self_energy_iw_im: {
    xLabel: 'i\u03c9<sub>n</sub>',
    yLabel: 'Im \u03A3(i\u03c9<sub>n</sub>)',
    axis: 'matsubara_freq'
  },
  greens_function_tau_re: {
    xLabel: '\u03c4<sub>n</sub>',
    yLabel: 'Re G(\u03c4<sub>n</sub>)',
    axis: 'tau'
  },
  greens_function_freq_im: {
    xLabel: '\u03c9',
    yLabel: 'Im G(\u03c9)',
    axis: 'frequencies'
  },
  self_energy_freq_im: {
    xLabel: '\u03c9',
    yLabel: 'Im \u03A3(\u03c9)',
    axis: 'frequencies'
  }
}

const useStyles = makeStyles(theme => ({
  plot: {height: '250px'}
}))

/**
 * Displays Greens functions data together with the provenance information.
 */
const GreensFunctions = React.memo(({
  data,
  provenance,
  className,
  'data-testid': testID,
  ...other
}) => {
  const theme = useTheme()
  const styles = useStyles()
  const [finalData, setFinalData] = useState(
    data || Object.fromEntries(Object.keys(types).map(type => [type, [data]]))
  )

  // Side effect that runs when the data that is displayed should change. By
  // running all this heavy stuff as a side effect, the first render containing
  // the placeholders etc. can be done as fast as possible.
  useEffect(() => {
    if (!data) {
      setFinalData(Object.fromEntries(Object.keys(types).map(type => [type, [data]])))
      return
    }

    const finalData = {}
    for (const [key, plotArray] of Object.entries(data)) {
      if (!plotArray) {
        finalData[key] = [plotArray]
        continue
      }
      const plots = plotArray.map(plot => {
        const {plotData, legendTitle, legendWidth} = plotDegeneracy(plot.x, plot.y, provenance)

        // Defining the plot
        const layout = {
          xaxis: {
            fixedrange: false,
            zeroline: true,
            title: {
              text: types[key].xLabel
            }
          },
          yaxis: {
            fixedrange: false,
            zeroline: true,
            title: {
              text: types[key].yLabel
            }
          },
          showlegend: true,
          legend: {
            x: 1,
            y: 1,
            title: {
              text: legendTitle,
              font: {
                size: 12
              }
            },
            xanchor: 'right',
            itemwidth: legendWidth
          }
        }
        return {
          type: key,
          plotData,
          layout,
          title: plot?.label || ''
        }
      })
      finalData[key] = plots
    }
    setFinalData(finalData)
  }, [data, provenance, theme])

  return <Box display="flex" flexDirection="column" height="100%" width="100%">
    {Object.entries(finalData).map(([key, plotArray]) => {
      if (!plotArray) plotArray = [plotArray]
      return plotArray.map((value, index) =>
        <Plot
          key={`${key}-${index}`}
          data={value === false ? false : value?.plotData}
          layout={value?.layout}
          floatTitle={value?.title}
          data-testid={`${testID}-${key}`}
          className={styles.plot}
          {...other}
        />
      )
    })}
    {(provenance && !Array.isArray(provenance)) && <Box flex="0 0 auto">
      <PropertyProvenanceList xs={12}>
        <PropertyProvenanceItem
          title="DMFT"
          data={provenance}
          path={"results.method.simulation.dmft"}
        />
      </PropertyProvenanceList>
    </Box>}
  </Box>
})

const plotType = PropTypes.oneOfType([PropTypes.arrayOf(PropTypes.object), PropTypes.oneOf([false, undefined])])
GreensFunctions.propTypes = {
  data: PropTypes.shape({
    greens_function_iw_im: plotType, // Real part of Green's functions in matsubara_freq
    self_energy_iw_im: plotType, // Imaginary part of self-energy in matsubara_freq
    greens_function_tau_re: plotType, // Real part of Green's functions in tau
    greens_function_freq_im: plotType, // Real part of Green's functions in frequencies
    self_energy_freq_im: plotType // Imaginary part of self-energy in frequencies
  }),
  provenance: PropTypes.object, // DMFT provenance
  className: PropTypes.string,
  'data-testid': PropTypes.string
}
GreensFunctions.defaultProps = {
  'data-testid': 'gfs'
}

/**
 * Given index and archive, resolves the data for plotting Green's functions.
 */
export function resolveGreensFunctions(properties, archive, pattern) {
  // Return undefined if property not available
  if (!properties.has('greens_functions_electronic')) return false

  // With the current metainfo for Green's functions, we cannot tell just with
  // the index data what type of plots will be available. Because of this, we
  // don't show anything until the archive is loaded.
  // TODO: This could be fixed by refactoring the greens functions data so that
  // each section would contain an indexed list of what plots it contains.
  if (!archive) return false

  // Resolve the final data once archive is here
  let gfReferences = archive?.results?.properties?.electronic?.greens_functions_electronic || []
  if (!Array.isArray(gfReferences)) gfReferences = [gfReferences]
  const greensFunctions = []
  for (const reference of gfReferences) {
    // Resolving Green's function quantities in tau
    const gfTau = {}
    const matchTau = reference.tau ? reference.tau.match(pattern) : undefined
    const pathTau = matchTau ? matchTau[2] : reference.tau
    const matchGreensFunctionTau = reference.greens_function_tau ? reference.greens_function_tau.match(pattern) : undefined
    const pathGreensFunctionTau = matchGreensFunctionTau ? matchGreensFunctionTau[2] : reference.greens_function_tau
    const sourceArchiveTau = matchTau
      ? (archive.m_ref_archives[matchTau[1]] || archive.m_ref_archives[reference.tau.split('#')[0]])
      : archive
    if (sourceArchiveTau && pathTau) {
      // @TODO ask why I have to flatten the axes (tau, matsubara_freq, frequencies) arrays
      gfTau.tau = resolveInternalRef(pathTau, sourceArchiveTau).flat()
      const greensFunctionTau = resolveInternalRef(pathGreensFunctionTau, sourceArchiveTau)
      gfTau.greens_function_tau_re = greensFunctionTau ? greensFunctionTau.re : undefined
    }

    // Resolving Green's function quantities in Matsubara frequencies
    const gfMatsFreq = {}
    const matchMatsFreq = reference.matsubara_freq ? reference.matsubara_freq.match(pattern) : undefined
    const pathMatsFreq = matchMatsFreq ? matchMatsFreq[2] : reference.matsubara_freq
    const matchGreensFunctionMatsFreq = reference.greens_function_iw ? reference.greens_function_iw.match(pattern) : undefined
    const pathGreensFunctionMatsFreq = matchGreensFunctionMatsFreq ? matchGreensFunctionMatsFreq[2] : reference.greens_function_iw
    const matchSelfEnergyMatsFreq = reference.self_energy_iw ? reference.self_energy_iw.match(pattern) : undefined
    const pathSelfEnergyMatsFreq = matchSelfEnergyMatsFreq ? matchSelfEnergyMatsFreq[2] : reference.self_energy_iw
    const sourceArchiveMatsFreq = matchMatsFreq
      ? (archive.m_ref_archives[matchMatsFreq[1]] || archive.m_ref_archives[reference.matsubara_freq.split('#')[0]])
      : archive
    if (sourceArchiveMatsFreq && pathMatsFreq) {
      gfMatsFreq.matsubara_freq = resolveInternalRef(pathMatsFreq, sourceArchiveMatsFreq).flat()
      const greensFunctionMatsFreq = resolveInternalRef(pathGreensFunctionMatsFreq, sourceArchiveMatsFreq)
      gfMatsFreq.greens_function_iw_im = greensFunctionMatsFreq ? greensFunctionMatsFreq.im : undefined
      const selfEnergyMatsFreq = resolveInternalRef(pathSelfEnergyMatsFreq, sourceArchiveMatsFreq)
      gfMatsFreq.self_energy_iw_im = selfEnergyMatsFreq ? selfEnergyMatsFreq.im : undefined
    }
    const gfImAxis = mergeObjects(gfTau, gfMatsFreq)

    // Resolving Green's function quantities in Matsubara frequencies
    const gfFreq = {}
    const matchFreq = reference.frequencies ? reference.frequencies.match(pattern) : undefined
    const pathFreq = matchFreq ? matchFreq[2] : reference.frequencies
    const matchGreensFunctionFreq = reference.greens_function_freq ? reference.greens_function_freq.match(pattern) : undefined
    const pathGreensFunctionFreq = matchGreensFunctionFreq ? matchGreensFunctionFreq[2] : reference.greens_function_freq
    const matchSelfEnergyFreq = reference.self_energy_freq ? reference.self_energy_freq.match(pattern) : undefined
    const pathSelfEnergyFreq = matchSelfEnergyFreq ? matchSelfEnergyFreq[2] : reference.self_energy_freq
    const sourceArchiveFreq = matchFreq
      ? (archive.m_ref_archives[matchFreq[1]] || archive.m_ref_archives[reference.frequencies.split('#')[0]])
      : archive
    if (sourceArchiveFreq && pathFreq) {
      gfFreq.frequencies = resolveInternalRef(pathFreq, sourceArchiveFreq).flat()
      const greensFunctionFreq = resolveInternalRef(pathGreensFunctionFreq, sourceArchiveFreq)
      gfFreq.greens_function_freq_im = greensFunctionFreq ? greensFunctionFreq.im : undefined
      const selfEnergyFreq = resolveInternalRef(pathSelfEnergyFreq, sourceArchiveFreq)
      gfFreq.self_energy_freq_im = selfEnergyFreq ? selfEnergyFreq.im : undefined
    }

    // Merging all resolved data into gf.
    const gf = mergeObjects(gfImAxis, gfFreq)
    gf.type = reference.type || undefined
    gf.label = reference.label || undefined
    gf.m_path = `${archive?.metadata?.entry_id}/data/results/properties/electronic/greens_functions_electronic`
    if (Object.keys(gf).length !== 0 && gf.type !== 'impurity') greensFunctions.push(gf)
  }

  // Store all data under one dictionary for plotting
  const data = {}
  for (const gf of greensFunctions) {
    for (const [type, value] of Object.entries(types)) {
      if (gf[type] && gf[value.axis]) {
        const oldData = data[type]
        if (!oldData) {
          data[type] = []
        }
        data[type].push({x: gf[value.axis], y: gf[type]})
      }
    }
  }
  return data
}

/**
 * Plots data in atom, orbital and spin space, taking into account degeneracies for the
 * labeling and titles.
 */
function plotDegeneracy(axisData, data, provenance) {
  const oneSpecies = data.length === 1
  const oneSpin = provenance.magnetic_state === 'paramagnetic'
  const atomLabels = alphabet
  let legendTitle = ''
  let legendWidth = 30
  const plotData = []
  for (let iAtom = 0; iAtom < data.length; ++iAtom) {
    const atomData = data[iAtom]
    if (oneSpin) { // We only show one spin data
      const spinData = atomData[0]
      for (let iOrb = 0; iOrb < spinData.length; ++iOrb) {
        plotData.push({
          x: axisData,
          y: spinData[iOrb],
          name: `${oneSpecies ? '' : `${atomLabels[iAtom]}`}d<sub>${iOrb}</sub>`,
          type: 'scatter'
        })
      }
      legendTitle = `${oneSpecies ? '<b>[orbital]</b>' : '<b>[atom][orbital]</b>'}`
      legendWidth = oneSpecies ? 30 : 40
    } else {
      for (let iSpin = 0; iSpin < atomData.length; ++iSpin) {
        const spinData = atomData[iSpin]
        for (let iOrb = 0; iOrb < spinData.length; ++iOrb) {
          plotData.push({
            x: axisData,
            y: spinData[iOrb],
            name: `${oneSpecies ? '' : `${atomLabels[iAtom]}`}d<sub>${iOrb}</sub>${oneSpin ? '' : iSpin ? '+' : '-'}`,
            type: 'scatter'
          })
        }
      }
      legendTitle = `${oneSpecies ? '<b>[orbital][spin]</b>' : '<b>[atom][orbital][spin]</b>'}`
      legendWidth = oneSpecies ? 40 : 60
    }
  }
  return {plotData, legendTitle, legendWidth}
}

export default withErrorHandler(gfError)(GreensFunctions)
