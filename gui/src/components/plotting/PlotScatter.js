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
import React, {useState, useEffect, useMemo, useCallback, forwardRef} from 'react'
import PropTypes from 'prop-types'
import { makeStyles, useTheme } from '@material-ui/core'
import { hasWebGLSupport } from '../../utils'
import * as d3 from 'd3'
import FilterTitle from '../search/FilterTitle'
import Plot from './Plot'
import { useHistory } from 'react-router-dom'
import { getUrl } from '../nav/Routes'

/**
 * A Plotly-based interactive scatter plot.
 */
const hasWebGL = hasWebGLSupport()
const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    height: '100%',
    boxSizing: 'border-box',
    display: 'grid',
    gridTemplateColumns: 'auto 1fr',
    gridTemplateRows: '1fr auto',
    marginBottom: theme.spacing(0.5)
  },
  plot: {
    gridColumn: 2,
    gridRow: 1,
    position: 'relative'
  },
  xaxis: {
    gridColumn: 2,
    gridRow: 2,
    display: 'flex',
    flexDirection: 'row',
    alignItems: 'center',
    justifyContent: 'center',
    height: '1rem',
    marginTop: theme.spacing(0.5)
  },
  yaxis: {
    gridColumn: 1,
    gridRow: 1,
    height: '100%',
    marginRight: theme.spacing(0.5),
    width: '1rem',
    display: 'flex',
    flexDirection: 'row',
    alignItems: 'center',
    justifyContent: 'center'
  },
  square: {
    gridColumn: 1,
    gridRow: 2
  },
  color: {
    gridColumn: 3,
    gridRow: 1,
    height: '100%',
    marginLeft: theme.spacing(0.5),
    width: '1rem',
    display: 'flex',
    flexDirection: 'row',
    alignItems: 'center',
    justifyContent: 'center'
  },
  colorlabel: {
    transform: 'rotate(90deg)'
  }
}))
const PlotScatter = React.memo(forwardRef((
{
  data,
  title,
  xAxis,
  yAxis,
  colorAxis,
  discrete,
  autorange,
  dragmode,
  onSelected,
  onDeselect,
  onNavigateToEntry,
  'data-testid': testID
}, canvas) => {
  const styles = useStyles()
  const theme = useTheme()
  const [finalData, setFinalData] = useState(!data ? data : undefined)
  const history = useHistory()

  // Side effect that runs when the data that is displayed should change. By
  // running all this heavy stuff within useEffect (instead of e.g. useMemo),
  // the first render containing the placeholders etc. can be done as fast as
  // possible.
  useEffect(() => {
    if (!data) {
      setFinalData(data)
      return
    }

    const hoverTemplate = (xLabel, yLabel, colorLabel, xUnit, yUnit, colorUnit) => {
      let template = `<b>Click to go to entry page</b>` +
        `<br>` +
        `${xLabel || ''}: %{x} ${xUnit === 'dimensionless' ? '' : xUnit}<br>` +
        `${yLabel || ''}: %{y} ${yUnit === 'dimensionless' ? '' : yUnit}<br>`
      if (colorLabel) {
        template = template +
          `${colorLabel || ''}: %{${discrete ? 'text' : 'text:.3'}} ${colorUnit === 'dimensionless' ? '' : colorUnit}<br>`
      }
      template = template + `<extra></extra>`
      return template
    }

    // If dealing with a quantized color, each group is separated into it's own
    // trace which has a legend as well.
    const traces = []
    if (colorAxis?.quantity && discrete) {
      const options = [...new Set(data.color)]
      const nOptions = options.length
      const scale = d3.scaleSequential([0, 1], d3.interpolateTurbo)
      const offset = 0.1
      for (const option of options) {
        const xArray = []
        const yArray = []
        const colorArray = []
        const entryIdArray = []
        for (let i = 0; i < data.color.length; ++i) {
          if (data.color[i] === option) {
            xArray.push(data.x[i])
            yArray.push(data.y[i])
            colorArray.push(data.color[i])
            entryIdArray.push(data.id[i])
          }
        }
        traces.push({
          x: xArray,
          y: yArray,
          entry_id: entryIdArray,
          name: option,
          text: colorArray,
          mode: 'markers',
          type: hasWebGL ? 'scattergl' : 'scatter',
          textposition: 'top center',
          showlegend: true,
          hovertemplate: hoverTemplate(
            xAxis.title,
            yAxis.title,
            colorAxis.title,
            xAxis.unit,
            yAxis.unit,
            ''
          ),
          marker: {
            size: 8,
            color: scale(offset + (1 - 2 * offset) * options.indexOf(option) / (nOptions - 1)),
            line: {
              color: theme.palette.grey[800],
              width: 1
            }
          }
        })
      }
    // When dealing with a continuous color, display a colormap
    } else if (colorAxis?.quantity && !discrete) {
      traces.push({
        x: data.x,
        y: data.y,
        color: data.color,
        text: data.color,
        entry_id: data.id,
        mode: 'markers',
        type: 'scattergl',
        textposition: 'top center',
        showlegend: false,
        hoverinfo: "text",
        hovertemplate: hoverTemplate(
          xAxis.title,
          yAxis.title,
          colorAxis.title,
          xAxis.unit,
          yAxis.unit,
          colorAxis.unit
        ),
        marker: {
          size: 8,
          color: data.color,
          colorscale: 'YlGnBu',
          line: {
            color: theme.palette.grey[800],
            width: 1
          },
          colorbar: {
            thickness: 20,
            ypad: 0,
            xpad: 5,
            tickfont: {
              family: theme.typography.fontFamily
            }
          }
        }
      })
    // When color is not set, all points are displayed in a single plot with
    // primary theme color.
    } else {
      traces.push({
        x: data.x,
        y: data.y,
        entry_id: data.id,
        mode: 'markers',
        type: 'scattergl',
        textposition: 'top center',
        showlegend: false,
        hoverinfo: "text",
        hovertemplate: hoverTemplate(
          xAxis.title,
          yAxis.title,
          '',
          xAxis.unit,
          yAxis.unit,
          ''
        ),
        marker: {
          size: 8,
          color: theme.palette.secondary.main,
          line: {
            color: theme.palette.grey[800],
            width: 1
          }
        }
      })
    }
    setFinalData(traces)
  }, [colorAxis?.quantity, colorAxis?.title, colorAxis?.unit, data, discrete, theme, xAxis.title, xAxis.unit, yAxis.title, yAxis.unit])

  const layout = useMemo(() => {
    return {
      dragmode: dragmode,
      hovermode: 'closest',
      hoverlabel: {
        bgcolor: theme.palette.grey[100],
        bordercolor: theme.palette.grey[100],
        font: {
          color: theme.palette.grey[800],
          family: theme.typography.fontFamily
        }
      },
      showlegend: true,
      legend: {
        x: 1,
        xanchor: 'right',
        y: 1
      },
      xaxis: {
        fixedrange: false,
        autorange: autorange
      },
      yaxis: {
        fixedrange: false,
        autorange: autorange
      },
      margin: {
        l: 8,
        r: 0,
        t: 8,
        b: 24
      }
    }
  // Any further changes to dragmode need to be handled through callbacks in
  // order to do the layout updates correctly. TODO: The plot component should
  // have a consistent way for handling layout changes: they should either
  // happen through a property or through some update function, but not through
  // both. This is a general problem in trying to 'reactify' a non-react library
  // like Plotly.
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [autorange])

  // Change dragmode
  useEffect(() => {
    canvas.current?.relayout && canvas.current.relayout((layout) => {
      return {
        ...layout,
        dragmode
      }
    })
  }, [dragmode, canvas])

  const handleClick = useCallback(d => {
    const pointIndex = d.points[0].pointIndex
    const entryId = d.points[0].data.entry_id[pointIndex]
    const path = `entry/id/${entryId}`
    onNavigateToEntry()
    history.push(getUrl(path))
  }, [history, onNavigateToEntry])

  return <div className={styles.root}>
    <div className={styles.yaxis}>
      <FilterTitle
        quantity={yAxis.quantity}
        label={yAxis.title}
        unit={yAxis.unit}
        variant="caption"
        rotation="up"
      />
    </div>
    <div className={styles.plot}>
      <Plot
        data={finalData}
        layout={layout}
        floatTitle={title}
        fixedMargins={true}
        autorange={autorange}
        onSelected={onSelected}
        onDeselect={onDeselect}
        disableDefaultActions
        throttleResize
        data-testid={testID}
        ref={canvas}
        onClick={handleClick}
      />
    </div>
    <div className={styles.square} />
    <div className={styles.xaxis}>
      <FilterTitle
        quantity={xAxis.quantity}
        label={xAxis.title}
        unit={xAxis.unit}
        variant="caption"
      />
    </div>
    {!discrete && colorAxis &&
      <div className={styles.color}>
        <FilterTitle
          rotation="down"
          quantity={colorAxis.quantity}
          unit={colorAxis.unit}
          label={colorAxis.title}
          description=""
          variant="caption"
        />
      </div>
    }
  </div>
}))

PlotScatter.propTypes = {
  data: PropTypes.object,
  title: PropTypes.string,
  xAxis: PropTypes.object, // Contains x-axis settings
  yAxis: PropTypes.object, // Contains y-axis settings
  colorAxis: PropTypes.object, // Contains colorbar settings
  discrete: PropTypes.bool,
  autorange: PropTypes.bool,
  dragmode: PropTypes.string,
  onSelected: PropTypes.func,
  onDeselect: PropTypes.func,
  onNavigateToEntry: PropTypes.func,
  'data-testid': PropTypes.string
}

PlotScatter.defaultProps = {
  unitX: 'dimensionless',
  unitY: 'dimensionless',
  unitColor: 'dimensionless'
}

export default PlotScatter
