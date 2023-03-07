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
import { useTheme, makeStyles } from '@material-ui/core/styles'
import { Box } from '@material-ui/core'
import Plot from '../plotting/Plot'
import {
    PropertyMethodologyItem,
    PropertyMethodologyList
} from '../entry/properties/PropertyCard'
import { withErrorHandler } from '../ErrorHandler'
import { mergeObjects, alphabet } from '../../utils'

/**
 * Displays Greens functions data as two plots (real + imag) together with the
 * methodology information.
 */
export const gfError = 'Could not load Greens functions.'
const types = ['regtau', 'imsiw']
const useStyles = makeStyles((theme) => ({
  regtau: {},
  imsiw: {}
}))
/*
 * Plots data in atom, orbital and spin space, taking into account degeneracies for the
 * labeling and titles.
 */
function plotDegeneracy(data, methodology, axisData) {
  const oneSpecies = data.length === 1
  const oneSpin = methodology.magnetic_state === 'paramagnetic'
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
const GreensFunctions = React.memo(({
  data,
  methodology,
  className,
  classes,
  'data-testid': testID,
  ...other
}) => {
  const theme = useTheme()
  const styles = useStyles({classes: classes})
  const [finalData, setFinalData] = useState(!data ? {regtau: data, imsiw: data} : data)
  const [finalLayout, setFinalLayout] = useState()

  // Side effect that runs when the data that is displayed should change. By
  // running all this heavy stuff as a side effect, the first render containing
  // the placeholders etc. can be done as fast as possible.
  useEffect(() => {
    if (!data) {
      setFinalData({regtau: data, imsiw: data})
      return
    }

    // Generate plot data for both regtau and imsiw
    const tracesAll = {}
    const layoutAll = {}
    for (const type of types) {
      let xLabel = ''
      let yLabel = ''
      let axisData = []
      let gfData = []
      if (type === 'regtau') {
        xLabel = '\u03c4<sub>n</sub>'
        yLabel = 'Re G(\u03c4<sub>n</sub>)'
        axisData = data.tau
        gfData = data.regtau
      } else {
        xLabel = 'i\u03c9<sub>n</sub>'
        yLabel = 'Im \u03A3(i\u03c9<sub>n</sub>)'
        axisData = data.iw
        gfData = data.imsiw
      }
      const initialLayout = {
        xaxis: {
          title: {
            text: xLabel
          },
          fixedrange: false,
          zeroline: true
        },
        yaxis: {
          title: {
            text: yLabel
          },
          fixedrange: false,
          zeroline: true
        }
      }

      // Defining the plot
      const plotting = plotDegeneracy(gfData, methodology, axisData)
      tracesAll[type] = plotting.plotData

      // Setting finalLayout: title for the legend
      const finalLayout = {
        showlegend: true,
        legend: {
          x: 1,
          y: 1,
          title: {
            text: plotting.legendTitle,
            font: {
              size: 12
            }
          },
          xanchor: 'right',
          itemwidth: plotting.legendWidth
        }
      }
      layoutAll[type] = mergeObjects(initialLayout, finalLayout)
    }
    setFinalLayout(layoutAll)
    setFinalData(tracesAll)
  }, [data, methodology, theme])

  return <Box display="flex" flexDirection="column" height="100%" width="100%">
    <Box flex="1 1 auto">
      {types.map((x) => <Plot
        key={x}
        data={finalData?.[x]}
        layout={finalLayout?.[x]}
        floatTitle={"Green's functions"}
        className={styles?.[x]}
        data-testid={`${testID}-${x}`}
        {...other}
      />)}
    </Box>
    {methodology && <Box flex="0 0 auto">
      <PropertyMethodologyList xs={12}>
        <PropertyMethodologyItem
          title="DMFT"
          data={methodology}
          path={"results.method.simulation.dmft"}
        />
      </PropertyMethodologyList>
    </Box>}
  </Box>
})
GreensFunctions.propTypes = {
  data: PropTypes.oneOfType([
    PropTypes.oneOf([false, undefined]), // false = no data, undefined = loading
    PropTypes.shape({
      tau: PropTypes.array, // Array of tau
      regtau: PropTypes.array, // Array of real_green_functions_tau
      iw: PropTypes.array, // Array of matsubara_freq
      imsiw: PropTypes.array // Array of imag_self_energy_iw
    })
  ]),
  methodology: PropTypes.object, // DMFT methodology
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}
GreensFunctions.defaultProps = {
  'data-testid': 'greens-functions'
}

export default withErrorHandler(gfError)(GreensFunctions)
