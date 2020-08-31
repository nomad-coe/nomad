import React, { useRef, useState, useEffect } from 'react'
import PropTypes from 'prop-types'
import { Select, MenuItem, Card, CardContent, CardHeader, makeStyles } from '@material-ui/core'
import * as d3 from 'd3'
import { scaleBand, scalePow } from 'd3-scale'
import { formatQuantity, nomadPrimaryColor, nomadSecondaryColor, nomadFontFamily } from '../config.js'

function split(array, cols) {
  if (cols === 1) {
    return [array]
  }
  const size = Math.ceil(array.length / cols)
  return [array.slice(0, size), ...split(array.slice(size), cols - 1)]
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
    textAlign: 'center',
    position: 'absolute',
    pointerEvents: 'none',
    opacity: 0
  },
  tooltipContent: {
    // copy of the material ui popper style
    display: 'inline-block',
    color: '#fff',
    padding: '4px 8px',
    fontSize: '0.625rem',
    lineHeight: '1.4em',
    borderRadius: '4px',
    backgroundColor: '#616161'
  }
}))
export default function Histogram({
  initialScale = 1, getValueLabel, numberOfValues, title, data, columns = 1, tooltips,
  onClick, selected, card
}) {
  onClick = onClick || (() => null)
  data = data || []
  numberOfValues = numberOfValues || data.length
  getValueLabel = getValueLabel || (value => value.name)
  title = title || 'Histogram'

  const classes = useStyles()
  const containerRef = useRef()
  const [scale, setScale] = useState(initialScale)

  const handleItemClicked = item => onClick(item)

  useEffect(() => {
    for (let i = data.length; i < numberOfValues; i++) {
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

    const x = scalePow().range([0, width]).exponent(scale)

    // we use at least the domain 0..1, because an empty domain causes a weird layout
    const max = d3.max(data, d => d.value) || 1
    x.domain([0, max])

    const rectColor = d => selected === d.key ? nomadPrimaryColor.dark : nomadSecondaryColor.light
    const textColor = d => selected === d.key ? '#FFF' : '#000'

    const container = d3.select(containerRef.current)
    const tooltip = container.select('.' + classes.tooltip)
      .style('width', width + 'px')
      .style('opacity', 0)
    const tooltipContent = container.select('.' + classes.tooltipContent)
    const svg = container.select('svg')
      .attr('width', containerWidth)
      .attr('height', height)

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
        .on('click', d => handleItemClicked(d))

      let item = items.enter()
        .append('g')
        .attr('class', 'item')
        .attr('display', d => getValueLabel(d) === '' ? 'none' : 'show')

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
        .style('cursor', 'pointer')
        .on('click', d => handleItemClicked(d))

      item
        .on('mouseover', function(d) {
          d3.select(this).select('.background')
            .style('opacity', 0.08)
          if (tooltips) {
            tooltip.transition()
              .duration(200)
              .style('opacity', 1)
            tooltip
              .style('left', i * (width + 12) + 'px')
              .style('top', (y(d.key) + 32) + 'px')
            tooltipContent.html(getValueLabel(d))
          }
        })
        .on('mouseout', function(d) {
          d3.select(this).select('.background')
            .style('opacity', 0)
          if (tooltips) {
            tooltip.transition()
              .duration(200)
              .style('opacity', 0)
          }
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
      <div className={classes.tooltipContent}></div>
    </div>
    <svg />
  </div>

  if (card) {
    return <Card classes={{root: classes.root}}>
      <CardHeader
        title={title}
        titleTypographyProps={{variant: 'body1'}}
        action={(
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
  card: PropTypes.bool
}
