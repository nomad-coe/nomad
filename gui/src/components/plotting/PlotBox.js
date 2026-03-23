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
import React, { useMemo, useCallback, forwardRef } from 'react'
import PropTypes from 'prop-types'
import { isNil } from 'lodash'
import * as d3 from 'd3'
import { useRecoilValue } from 'recoil'
import { makeStyles, useTheme } from '@material-ui/core'
import { useHistory } from 'react-router-dom'
import { getUrl } from '../nav/Routes'
import { DefinitionTitle } from '../DefinitionTitle'
import Placeholder from '../visualization/Placeholder'
import { guiState } from '../GUIMenu'
import Plot from './Plot'
import { getScatterPlotHoverTemplate, getAxisType } from './common'

/**
 * A Plotly-based box plot. Renders pre-computed box plot statistics
 * (min, q1, median, q3, max) for one or more groups.
 */
const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    height: '100%',
    boxSizing: 'border-box',
    display: 'grid',
    marginBottom: theme.spacing(0.5)
  },
  plot: {
    position: 'relative'
  },
  xaxis: {
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
  }
}))

const titleClasses = {noWrap: true, maxWidth: '20em'}
const boxTraceConfig = {
  type: 'box',
  boxpoints: false,
  hoveron: 'boxes',
  hoverinfo: 'y',
  showlegend: false
}

// Offset of scatter points to the left of the box center and the half-width
// of the jitter band.
const pointOffset = -0.35
const jitterWidth = 0.02

/**
 * Simple deterministic hash-based jitter so the point positions are stable
 * across re-renders for the same data.
 */
function jitter(index) {
  let h = index * 2654435761 >>> 0 // Knuth multiplicative hash
  h = ((h >> 16) ^ h) * 0x45d9f3b >>> 0
  return pointOffset + (((h % 1000) / 999) - 0.5) * 2 * jitterWidth
}

const PlotBox = React.memo(forwardRef(({
  data,
  title,
  xAxis,
  yAxis,
  colorAxis,
  discrete,
  show_points,
  autorange,
  dragMode,
  onNavigateToEntry,
  'data-testid': testID
}, canvas) => {
  const hasYLabel = !!yAxis?.title
  const hasXLabel = !!xAxis?.title
  const styles = useStyles()
  const theme = useTheme()
  const history = useHistory()
  const aggIndicator = useRecoilValue(guiState('aggIndicator'))
  const hasLabels = data && data.some(item => item.label != null)
  const hasSubgroups = data && data.some(item => item.subgroups && item.subgroups.length > 0)
  const hasDiscreteColor = discrete && colorAxis?.search_quantity
  const hasContinuousColor = !discrete && colorAxis?.search_quantity

  // Stable mapping from labels to numeric x positions. When subgroups exist, every
  // (group, subgroup) pair gets its own sequential position so that all boxes are placed
  // side-by-side. Group tick labels are centered across the range of positions belonging
  // to that group.
  const { labelToPos, subPosMap, tickVals, tickText } = useMemo(() => {
    const labelToPos = {}
    if (!data) return { labelToPos, subPosMap: null, tickVals: [], tickText: [] }

    if (hasSubgroups) {
      const subPosMap = new Map()
      let pos = 0
      const groupRanges = []
      for (const group of data) {
        const groupKey = group.label ?? '__default__'
        const startPos = pos
        for (const sub of (group.subgroups || [])) {
          const key = `${groupKey}::${sub.label}`
          if (!subPosMap.has(key)) {
            subPosMap.set(key, pos)
            pos++
          }
        }
        if (pos > startPos) {
          groupRanges.push({ label: group.label, startPos, endPos: pos - 1 })
        }
      }
      return {
        labelToPos,
        subPosMap,
        tickVals: groupRanges.map(r => (r.startPos + r.endPos) / 2),
        tickText: groupRanges.map(r => r.label ?? '')
      }
    }

    // Non-subgroup case: one position per group label
    const labelOrder = []
    if (hasLabels) {
      for (const group of data) {
        if (group.label != null && !(group.label in labelToPos)) {
          labelToPos[group.label] = labelOrder.length
          labelOrder.push(group.label)
        }
      }
    }
    return {
      labelToPos,
      subPosMap: null,
      tickVals: labelOrder.map((_, i) => i),
      tickText: labelOrder
    }
  }, [data, hasLabels, hasSubgroups])

  // Build Plotly traces from the group objects. Each group has raw y/color/id arrays. Box
  // statistics are computed by Plotly from raw values. Individual data points are
  // rendered as a separate scatter trace so that hover templates and color encoding
  // (colorbar / legend) can be controlled independently.
  const traces = useMemo(() => {
    if (!data || data.length === 0) return []

    // Subgroup helpers: compute unique subgroup labels, a color scale, and an offset
    // function that places each subgroup's box side-by-side within a group, centered on
    // the group x position.
    let subLabels, nSub, subColorScale
    if (hasSubgroups) {
      subLabels = [...new Set(data.flatMap(d => (d.subgroups || []).map(s => s.label)).filter(l => l != null))]
      nSub = subLabels.length
      const scale = d3.scaleSequential([0, 1], d3.interpolateTurbo)
      const cOffset = 0.1
      subColorScale = (subIdx) => scale(cOffset + (1 - 2 * cOffset) * subIdx / Math.max(nSub - 1, 1))
    }

    const hasColor = hasSubgroups
      ? data.some(g => (g.subgroups || []).some(s => s.color && s.color.length > 0))
      : data.some(g => g.color && g.color.length > 0)

    // Box traces
    const boxTraces = []
    if (hasSubgroups) {
      const seenLegend = new Set()
      for (const group of data) {
        const groupKey = group.label ?? '__default__'
        for (const sub of (group.subgroups || [])) {
          const key = `${groupKey}::${sub.label}`
          const xPos = subPosMap.get(key) ?? 0
          const subIdx = subLabels.indexOf(sub.label)
          const color = subColorScale(subIdx)
          const showLegend = !seenLegend.has(sub.label)
          seenLegend.add(sub.label)
          boxTraces.push({
            ...boxTraceConfig,
            x: sub.y.map(() => xPos),
            y: sub.y,
            name: sub.label || '',
            legendgroup: sub.label || '',
            showlegend: showLegend,
            marker: { color },
            line: { color }
          })
        }
      }
    } else if (hasLabels) {
      for (const group of data) {
        const pos = labelToPos[group.label] ?? 0
        boxTraces.push({
          ...boxTraceConfig,
          x: group.y.map(() => pos),
          y: group.y,
          name: group.label || '',
          marker: { color: theme.palette.primary.light },
          line: { color: theme.palette.primary.light }
        })
      }
    } else {
      const allY = data.flatMap(g => g.y)
      boxTraces.push({
        ...boxTraceConfig,
        x: allY.map(() => 0),
        y: allY,
        marker: { color: theme.palette.primary.light },
        line: { color: theme.palette.primary.light }
      })
    }

    // Scatter traces
    const scatterTraces = []
    if (show_points) {
      const flatX = []
      const flatY = []
      const flatColor = []
      const flatId = []
      const flatSubLabel = []
      let globalIdx = 0
      const hoverTemplate = getScatterPlotHoverTemplate(undefined, yAxis?.title, colorAxis?.title, discrete)

      for (const group of data) {
        if (hasSubgroups) {
          if (!group.subgroups) continue
          const groupKey = group.label ?? '__default__'
          for (const sub of group.subgroups) {
            const key = `${groupKey}::${sub.label}`
            const xPos = subPosMap.get(key) ?? 0
            for (let i = 0; i < sub.y.length; i++) {
              flatX.push(xPos + jitter(globalIdx))
              flatY.push(sub.y[i])
              if (sub.color) flatColor.push(sub.color[i])
              flatId.push(sub.id[i])
              flatSubLabel.push(sub.label)
              globalIdx++
            }
          }
        } else {
          const pos = hasLabels ? (labelToPos[group.label] ?? 0) : 0
          for (let i = 0; i < group.y.length; i++) {
            flatX.push(pos + jitter(globalIdx))
            flatY.push(group.y[i])
            if (group.color) flatColor.push(group.color[i])
            flatId.push(group.id[i])
            globalIdx++
          }
        }
      }

      if (hasColor && hasDiscreteColor) {
        // Discrete color: one scatter trace per category with legend
        const options = [...new Set(flatColor)]
        const nOptions = options.length
        const scale = d3.scaleSequential([0, 1], d3.interpolateTurbo)
        const offset = 0.1
        for (const option of options) {
          const xArr = []
          const yArr = []
          const colorArr = []
          const idArr = []
          for (let i = 0; i < flatColor.length; i++) {
            if (flatColor[i] === option) {
              xArr.push(flatX[i])
              yArr.push(flatY[i])
              colorArr.push(flatColor[i])
              idArr.push(flatId[i])
            }
          }
          if (yArr.length === 0) continue
          scatterTraces.push({
            x: xArr,
            y: yArr,
            entry_id: idArr,
            name: option,
            text: colorArr,
            mode: 'markers',
            type: 'scatter',
            showlegend: true,
            legendgroup: option,
            hovertemplate: hoverTemplate,
            marker: {
              size: 4,
              color: scale(offset + (1 - 2 * offset) * options.indexOf(option) / (nOptions - 1)),
              line: { color: theme.palette.grey[800], width: 0.5 }
            }
          })
        }
      } else if (hasColor && hasContinuousColor) {
        // Continuous color: single scatter trace with colorbar
        scatterTraces.push({
          x: flatX,
          y: flatY,
          entry_id: flatId,
          text: flatColor,
          mode: 'markers',
          type: 'scatter',
          showlegend: false,
          hovertemplate: hoverTemplate,
          marker: {
            size: 4,
            color: flatColor,
            colorscale: 'YlGnBu',
            line: { color: theme.palette.grey[800], width: 0.5 },
            colorbar: {
              thickness: 20,
              ypad: 0,
              xpad: 5,
              tickfont: { family: theme.typography.fontFamily }
            }
          }
        })
      } else if (hasSubgroups) {
        // No color axis with subgroups: one trace per subgroup label
        subLabels.forEach((subLabel, s) => {
          const subX = []
          const subY = []
          const subId = []
          for (let i = 0; i < flatSubLabel.length; i++) {
            if (flatSubLabel[i] === subLabel) {
              subX.push(flatX[i])
              subY.push(flatY[i])
              subId.push(flatId[i])
            }
          }
          if (subY.length === 0) return
          scatterTraces.push({
            x: subX,
            y: subY,
            entry_id: subId,
            mode: 'markers',
            type: 'scatter',
            name: subLabel,
            showlegend: false,
            legendgroup: subLabel,
            hovertemplate: hoverTemplate,
            marker: {
              size: 4,
              color: subColorScale(s),
              line: { color: theme.palette.grey[800], width: 0.5 }
            }
          })
        })
      } else {
        // No color axis, no subgroups: single trace
        scatterTraces.push({
          x: flatX,
          y: flatY,
          entry_id: flatId,
          mode: 'markers',
          type: 'scatter',
          showlegend: false,
          hovertemplate: hoverTemplate,
          marker: {
            size: 4,
            color: theme.palette.primary.main,
            line: { color: theme.palette.grey[800], width: 0.5 }
          }
        })
      }
    }

    return [...boxTraces, ...scatterTraces]
  }, [data, hasLabels, hasSubgroups, hasDiscreteColor, hasContinuousColor, show_points, theme, colorAxis, discrete, yAxis?.title, labelToPos, subPosMap])

  // Click handler for navigating to entry page
  const handleClick = useCallback(d => {
    const point = d.points?.[0]
    if (!point) return
    const entryId = point.data.entry_id?.[point.pointIndex]
    if (!entryId) return
    const path = `entry/id/${entryId}`
    onNavigateToEntry?.()
    history.push(getUrl(path))
  }, [history, onNavigateToEntry])

  const layout = useMemo(() => ({
    hovermode: 'closest',
    dragmode: dragMode || 'zoom',
    hoverlabel: {
      bgcolor: theme.palette.grey[100],
      bordercolor: theme.palette.grey[100],
      font: {
        color: theme.palette.grey[800],
        family: theme.typography.fontFamily
      }
    },
    showlegend: !!hasDiscreteColor || !!hasSubgroups,
    legend: {
      x: 1,
      xanchor: 'right',
      y: 1
    },
    boxmode: hasSubgroups ? 'overlay' : undefined,
    xaxis: {
      autorange: true,
      fixedrange: false,
      ...(hasLabels || hasSubgroups
        ? {
          tickvals: tickVals,
          ticktext: tickText
        }
        : {showticklabels: false})
    },
    yaxis: {
      type: getAxisType(yAxis.dtype, yAxis.scale),
      fixedrange: false,
      autorange: autorange !== false,
      zeroline: false
    },
    margin: {
      l: 12,
      r: 0,
      t: 8,
      b: 24
    }
  }), [autorange, yAxis?.dtype, yAxis?.scale, dragMode, theme, hasLabels, hasDiscreteColor, hasSubgroups, tickVals, tickText])

  // Determine the plot content based on data state
  let plotContent
  if (isNil(data) && aggIndicator === 'on') {
    plotContent = <Placeholder
      variant="rect"
      data-testid={`${testID}-placeholder`}
      margin={0}
    />
  } else {
    plotContent = <Plot
      data={traces}
      layout={layout}
      floatTitle={title || undefined}
      fixedMargins={false}
      autorange={autorange}
      disableDefaultActions
      throttleResize
      data-testid={testID}
      ref={canvas}
      onClick={handleClick}
    />
  }

  return <div
    className={styles.root}
    style={{
      gridTemplateColumns: hasYLabel ? 'auto 1fr' : '1fr',
      gridTemplateRows: hasXLabel ? '1fr auto' : '1fr'
    }}
  >
    {hasYLabel && <div className={styles.yaxis}>
      <DefinitionTitle
        label={yAxis.title}
        description={yAxis?.description}
        variant="subtitle2"
        rotation="up"
        classes={titleClasses}
      />
    </div>}
    <div className={styles.plot}>
      {plotContent}
    </div>
    {hasXLabel && <>
      {hasYLabel && <div className={styles.square} />}
      <div className={styles.xaxis}>
        <DefinitionTitle
          label={xAxis.title}
          description={xAxis?.description}
          variant="subtitle2"
          classes={titleClasses}
        />
      </div>
    </>}
    {show_points && hasContinuousColor && colorAxis &&
      <div className={styles.color}>
        <DefinitionTitle
          label={colorAxis.title}
          description={colorAxis.description}
          rotation="down"
          variant="subtitle2"
          classes={titleClasses}
        />
      </div>
    }
  </div>
}))

PlotBox.propTypes = {
  // An array of group objects. Each group can optionally contain a list of subgroups for
  // hierarchical data.
  data: PropTypes.arrayOf(PropTypes.shape({
    label: PropTypes.string,
    y: PropTypes.array,
    color: PropTypes.array,
    id: PropTypes.array,
    subgroups: PropTypes.arrayOf(PropTypes.shape({
      label: PropTypes.string,
      y: PropTypes.array,
      color: PropTypes.array,
      id: PropTypes.array
    }))
  })),
  title: PropTypes.string,
  xAxis: PropTypes.object,
  yAxis: PropTypes.object,
  colorAxis: PropTypes.object,
  discrete: PropTypes.bool,
  show_points: PropTypes.bool,
  autorange: PropTypes.bool,
  dragMode: PropTypes.string,
  onNavigateToEntry: PropTypes.func,
  'data-testid': PropTypes.string
}

export default PlotBox
