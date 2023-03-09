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
import React, { useRef, useState, useEffect } from 'react'
import PropTypes from 'prop-types'
import { Select, MenuItem, Card, CardContent, CardHeader, makeStyles, Tooltip } from '@material-ui/core'
import * as d3 from 'd3'
import { scaleBand, scalePow } from 'd3-scale'
import { formatQuantity, nomadPrimaryColor, nomadSecondaryColor, nomadFontFamily, nomadTheme } from '../config.js'

function split(array, cols) {
  if (cols === 1) {
    return [array]
  }
  const size = Math.ceil(array.length / cols)
  return [array.slice(0, size), ...split(array.slice(size), cols - 1)]
}

function isSelected(d, selected, multiple) {
  // Determine if the value has been selected
  let isSelected
  if (multiple) {
    if (selected === undefined) {
      isSelected = false
    } else {
      if (Array.isArray(selected)) {
        const selections = new Set(selected)
        isSelected = selections.has(d.key)
      } else {
        isSelected = selected === d.key
      }
    }
  } else {
    isSelected = selected === d.key
  }
  return isSelected
}

const useStyles = makeStyles(theme => ({
  root: {
    marginTop: theme.spacing(2)
  },
  content: {
    paddingTop: 0,
    position: 'relative'
  },
  tooltip: {
    width: '100%',
    height: '100%'
  },
  tooltipContent: {
    position: 'absolute',
    fontSize: nomadTheme.overrides.MuiTooltip.tooltip.fontSize,
    fontWeight: nomadTheme.overrides.MuiTooltip.tooltip.fontWeight,
    // backgroundColor: '#ffbb00' // Uncomment for debugging tooltips
    display: 'none',
    zIndex: 1,
    cursor: 'pointer'
  },
  canvas: {
    zIndex: 2
  }
}))
export default function Histogram({
  initialScale = 1, getValueLabel, numberOfValues, title, data, columns = 1, tooltips,
  onClick, selected, card, multiple
}) {
  onClick = onClick || (() => null)
  data = data || []
  numberOfValues = numberOfValues || data.length
  getValueLabel = getValueLabel || (value => value.name)
  title = title || 'Histogram'

  const [tooltipHTML, setTooltipHTML] = useState('')
  const [item, setItem] = useState()
  const [scale, setScale] = useState(initialScale)
  const classes = useStyles()
  const containerRef = useRef()

  const handleItemClicked = item => {
    if (item.value !== 0) {
      onClick(item)
    }
  }

  useEffect(() => {
    // TODO add proper treatment of not processed on server side and processing
    const numberOfValuesToRender = numberOfValues
    // if (data[data.length - 1] && data[data.length - 1].key === 'not processed') {
    //   if (data[data.length - 2].key === 'unavailable') {
    //     data[data.length - 2].value += data[data.length - 1].value
    //     data.pop()
    //     numberOfValuesToRender -= 1
    //   }
    // }

    for (let i = data.length; i < numberOfValuesToRender; i++) {
      data.push({key: `empty${i}`, name: '', value: 0})
    }

    const columnSize = Math.ceil(data.length / columns)
    for (let i = data.length; i < columnSize * columns; i++) {
      data.push({key: `empty${i}`, name: '', value: 0})
    }
    const columnsData = split(data, columns)

    const containerWidth = containerRef.current.offsetWidth
    const width = containerWidth / columns - (12 * (columns - 1))
    const height = columnSize * 32
    const padding = 16

    const x = scalePow().range([0, width]).exponent(scale)

    // we use at least the domain 0..1, because an empty domain causes a weird layout
    const max = d3.max(data, d => d.value) || 1
    x.domain([0, max])

    const rectColor = d => isSelected(d, selected, multiple) ? nomadSecondaryColor.main : nomadPrimaryColor.light
    const textColor = d => {
      if (d.value === 0) {
        return '#999'
      }
      return isSelected(d, selected, multiple) ? '#FFF' : '#000'
    }

    const container = d3.select(containerRef.current)
    const tooltipContent = container.select('.' + classes.tooltipContent)
    const svg = container.select('svg')
      .attr('width', containerWidth)
      .attr('height', height)
    tooltipContent
      .style('width', width + 'px')
      .style('height', 27 + 'px')

    const columnsG = svg
      .selectAll('.column')
      .data(columnsData.map((_, i) => `column${i}`))

    columnsG.exit().remove()
    columnsG
      .enter()
      .append('g')
      .attr('id', d => d)
      .attr('class', 'column')
      .attr('transform', (d, i) => `translate(${i * (width + 12)}, 0)`)

    columnsData.forEach((data, i) => {
      const y = scaleBand().rangeRound([0, height]).padding(0.1)
      y.domain(data.map(d => d.key))

      const items = svg.select('#column' + i)
        .selectAll('.item')
        .data(data, d => d.key)

      items.exit().remove()

      items
        .on('click', handleItemClicked)

      let item = items.enter()
        .append('g')
        .attr('class', 'item')
        .attr('display', d => getValueLabel(d) === '' ? 'none' : 'show')
        .style('cursor', 'pointer')
        .on('click', handleItemClicked)

      item
        .append('rect')
        .attr('x', x(0))
        .attr('y', d => y(d.key))
        .attr('width', width)
        .attr('class', 'background')
        .style('opacity', 0)
        .attr('height', y.bandwidth())

      item
        .append('rect')
        .attr('class', 'bar')
        .attr('x', x(0))
        .attr('y', d => y(d.key))
        .attr('width', d => x(d.value) - x(0))
        .attr('height', y.bandwidth())
        .style('fill', rectColor)
        // .style('stroke', '#000')
        // .style('stroke-width', '1px')
        .style('shape-rendering', 'geometricPrecision')

      item
        .append('text')
        .attr('class', 'name')
        .attr('dy', '.75em')
        .attr('x', x(0) + 4)
        .attr('y', d => y(d.key) + 4)
        .attr('text-anchor', 'start')
        .style('fill', textColor)
        .style('font-family', nomadFontFamily)
        .text(d => getValueLabel(d))

      item
        .append('text')
        .attr('class', 'value')
        .attr('dy', y.bandwidth())
        .attr('y', d => y(d.key) - 4)
        .attr('x', d => width - 4)
        .attr('text-anchor', 'end')
        .style('fill', textColor)
        .style('font-family', nomadFontFamily)
        .text(d => formatQuantity(d.value))

      item
        .on('mouseenter', function(d) {
          setItem(d)
          if (tooltips) {
            if (d.tooltip) {
              tooltipContent
                .style('left', i * (width) + padding + 'px')
                .style('top', (y(d.key)) + 'px')
                .style('display', 'block')
              setTooltipHTML(d.tooltip)
            }
          }
        })
        .on('mouseleave', function(d) {
        })

      item = items.transition(d3.transition().duration(500))

      item
        .select('.bar')
        .attr('y', d => y(d.key))
        .attr('width', d => x(d.value) - x(0))
        .attr('height', y.bandwidth())
        .style('fill', rectColor)

      item
        .select('.name')
        .text(d => getValueLabel(d))
        .attr('y', d => y(d.key) + 4)
        .style('fill', textColor)

      item
        .select('.value')
        .text(d => formatQuantity(d.value))
        .attr('y', d => y(d.key) - 4)
        .attr('x', width - 4)
        .style('fill', textColor)
    })
  })

  const chart = <div ref={containerRef}>
    <div className={classes.tooltip}>
      <Tooltip
        interactive
        placement='right'
        title={tooltipHTML}
      >
        <div className={classes.tooltipContent} onClick={() => { handleItemClicked(item) }}>
        </div>
      </Tooltip>
      <svg className={classes.canvas}></svg>
    </div>
  </div>

  if (card) {
    return <Card classes={{root: classes.root}}>
      <CardHeader
        title={title}
        titleTypographyProps={{variant: 'body1'}}
        action={(
          <Tooltip title="Select the power of the scale">
            <Select
              value={scale}
              onChange={(event) => setScale(event.target.value)}
              displayEmpty
              name="scale power"
            >
              <MenuItem value={1}>linear</MenuItem>
              <MenuItem value={0.5}>1/2</MenuItem>
              <MenuItem value={0.25}>1/4</MenuItem>
              <MenuItem value={0.125}>1/8</MenuItem>
            </Select>
          </Tooltip>
        )}
      />
      <CardContent classes={{root: classes.content}}>
        {chart}
      </CardContent>
    </Card>
  } else {
    return chart
  }
}
Histogram.propTypes = {
  /**
   * An optional title for the chart. If no title is given, the Histogram is used.
   */
  title: PropTypes.string,
  /**
   * An optional scale power that is used as the initial scale before the user
   * changes it. Default is 1 (linear scale).
   */
  initialScale: PropTypes.number,
  /**
   * The data. Must be an array of object with keys key (unique key), value (number),
   * and name. Name is only necessary if no getValueLabel is given.
   */
  data: PropTypes.arrayOf(PropTypes.object),
  /**
   * The maximum number of values. This is used to fix the histograms size. Otherwise,
   * the size is determined by the required space to render the existing values.
   */
  numberOfValues: PropTypes.number,
  /**
   * The number of columns that the values should be rendered in. Default is one.
   */
  columns: PropTypes.number,
  /**
   * Set to true to enable tooltips for each value.
   */
  tooltips: PropTypes.bool,
  /**
   * This callback is called if a value is clicked.
   */
  onClick: PropTypes.func,
  /**
   * The given value will be highlighted as selected.
   */
  selected: PropTypes.any,
  /**
   * A function that determined the label of a value.
   */
  getValueLabel: PropTypes.func,
  /**
   * If true the chart is wrapped in a MUI paper card with title and form.
   */
  card: PropTypes.bool,
  /**
   * If true, multiple values may be selected.
   */
  multiple: PropTypes.bool
}

Histogram.defaultProps = {
  multiple: false
}
