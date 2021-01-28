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
import React, { useState, useCallback, useMemo } from 'react'
import PropTypes from 'prop-types'
import {
  Box,
  Typography,
  useTheme
} from '@material-ui/core'
import { RecoilRoot } from 'recoil'
import { makeStyles } from '@material-ui/core/styles'
import Plot from '../visualization/Plot'
import Structure from '../visualization/Structure'
import { ErrorHandler } from '../ErrorHandler'

function GeoOptOverview({data, className, classes}) {
  const [step, setStep] = useState(0)

  // Styles
  const useStyles = makeStyles((theme) => {
    return {
      root: {
        display: 'flex',
        width: '100%'
      },
      energies: {
        flex: '1 1 66.6%'
      },
      structure: {
        flex: '1 1 33.3%'
      }
    }
  })
  const style = useStyles(classes)
  const theme = useTheme()
  const plotData = useMemo(() => {
    return [{
      x: [...Array(data.energies.length).keys()],
      y: data.energies,
      type: 'scatter',
      line: {
        color: theme.palette.primary.main,
        width: 2
      }
    }]
  }, [data, theme])

  const plotLayout = useMemo(() => {
    return {
      hovermode: 'x',
      hoverdistance: 100,
      spikedistance: 1000,
      xaxis: {
        title: 'Step number',
        tickmode: 'auto',
        autorange: true,
        zeroline: false,
        showspikes: true,
        spikethickness: 2,
        spikedash: 'dot',
        spikecolor: '#999999',
        spikemode: 'across' },
      yaxis: {
        title: 'Energy change (eV)',
        autorange: true,
        zeroline: false
      }
    }
  }, [])

  // Handles hover event on the plot to update the currently shown structure
  const handleHover = useCallback((event) => {
    setStep(event.points[0].x)
  }, [])

  return (
    <RecoilRoot>
      <Box className={style.root}>
        <Box className={style.energies}>
          <Typography variant="subtitle1" align='center'>Energy convergence</Typography>
          <ErrorHandler message='Could not load energies.'>
            <Plot
              data={plotData}
              layout={plotLayout}
              aspectRatio={1.5}
              onHover={handleHover}
              floatTitle="Energy convergence"
            >
            </Plot>
          </ErrorHandler>
        </Box>
        <Box className={style.structure}>
          <Typography variant="subtitle1" align='center'>Optimization trajectory</Typography>
          <ErrorHandler message='Could not load structure.'>
            <Structure
              system={data.structures[step]}
              aspectRatio={0.75}
              options={{view: {fitMargin: 0.75}}}
              positionsOnly={true}
            ></Structure>
          </ErrorHandler>
        </Box>
      </Box>
    </RecoilRoot>
  )
}

GeoOptOverview.propTypes = {
  data: PropTypes.object,
  className: PropTypes.string,
  classes: PropTypes.object
}

export default GeoOptOverview
