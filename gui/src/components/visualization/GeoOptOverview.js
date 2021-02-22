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
import React, { useCallback, useMemo } from 'react'
import { Subject } from 'rxjs'
import PropTypes from 'prop-types'
import {
  Box,
  Typography,
  useTheme
} from '@material-ui/core'
import { makeStyles } from '@material-ui/core/styles'
import Plot from '../visualization/Plot'
import { Structure } from '../visualization/Structure'
import { ErrorHandler, withErrorHandler } from '../ErrorHandler'

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

function GeoOptOverview({data, className, classes}) {
  // RxJS subject for efficiently propagating changes in structure information
  const positionsSubject = useMemo(() => new Subject(), [])

  // Styles
  const style = useStyles(classes)
  const theme = useTheme()

  const plotData = useMemo(() => {
    let steps = [...Array(data.structures.length).keys()]
    let energies = data.energies
    const diff = [0]
    for (let i = 0; i < energies.length - 1; ++i) {
      diff.push(energies[i + 1] - energies[i])
    }
    const traces = [
      {
        x: steps,
        y: diff,
        hoverinfo: 'x',
        name: 'Change per step',
        type: 'scatter',
        line: {
          color: theme.palette.secondary.main,
          width: 2
        }
      },
      {
        x: steps,
        y: energies,
        hoverinfo: 'x',
        name: 'Total change',
        type: 'scatter',
        line: {
          color: theme.palette.primary.main,
          width: 2
        }
      }
    ]
    if (data?.energy_change_criteria) {
      traces.push({
        x: [steps[0], steps[steps.length - 1]],
        y: [data?.energy_change_criteria, data?.energy_change_criteria],
        hoverinfo: 'x',
        name: 'Convergence criteria',
        type: 'line',
        line: {
          color: '#999',
          width: 1,
          dash: '10px,10px'
        }
      })
    }
    return traces
  }, [data, theme])

  const plotLayout = useMemo(() => {
    return {
      hovermode: 'x',
      hoverdistance: 100,
      spikedistance: 1000,
      showlegend: true,
      legend: {
        x: 0,
        xanchor: 'left',
        y: 1,
        bgcolor: 'rgba(255, 255, 255, 0.9)'
      },
      xaxis: {
        showexponent: 'first',
        title: {
          text: 'Step number'
        },
        tickmode: 'auto',
        tickformat: ',d',
        autorange: false,
        range: [0, data.structures.length - 1],
        zeroline: false,
        showspikes: true,
        spikethickness: 2,
        spikedash: 'dot',
        spikecolor: '#999999',
        spikemode: 'across' },
      yaxis: {
        title: {
          text: 'Energy (eV)'
        },
        autorange: true,
        zeroline: false
      }
    }
  }, [data])

  // Handles hover event on the plot to update the currently shown structure
  const handleHover = useCallback((event) => {
    positionsSubject.next(data.structures[event.points[0].x].positions)
  }, [data, positionsSubject])

  return (
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
        <Structure
          system={data.structures[0]}
          aspectRatio={0.75}
          options={{view: {fitMargin: 0.75}}}
          positionsSubject={positionsSubject}
        ></Structure>
      </Box>
    </Box>
  )
}

GeoOptOverview.propTypes = {
  data: PropTypes.object,
  className: PropTypes.string,
  classes: PropTypes.object
}

export default withErrorHandler(GeoOptOverview, 'Could not load geometry optimization data.')
