import React, { useContext, useState, useEffect, useRef, useLayoutEffect, useCallback } from 'react'
import PropTypes from 'prop-types'
import { Select, MenuItem, Card, CardHeader, CardContent, makeStyles } from '@material-ui/core'
import Grid from '@material-ui/core/Grid'
import TextField from '@material-ui/core/TextField'
import * as d3 from 'd3'
import { scaleTime, scalePow } from 'd3-scale'
import { nomadSecondaryColor } from '../../config.js'
import { searchContext, Dates } from './SearchContext'

const useStyles = makeStyles(theme => ({
  root: {
    marginTop: theme.spacing(2)
  },
  header: {
    paddingBottom: 0
  },
  content: {
    paddingTop: 0,
    position: 'relative',
    height: 250
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
    fontFamily: '"Roboto", "Helvetica", "Arial", sans-serif',
    lineHeight: '1.4em',
    borderRadius: '4px',
    backgroundColor: '#616161'
  }
}))

export default function UploadsHistogram({title = 'Uploads over time', initialScale = 1, tooltips}) {
  const classes = useStyles()
  const containerRef = useRef()
  const fromTimeFieldRef = useRef()
  const untilTimeFieldRef = useRef()
  const [scale, setScale] = useState(initialScale)
  const {response, query, setQuery, domain, setDateHistogram} = useContext(searchContext)

  useEffect(() => {
    setDateHistogram(true)
    return () => {
      setDateHistogram(false)
    }
  }, [setDateHistogram])

  useLayoutEffect(() => {
    fromTimeFieldRef.current.value = Dates.FormDate(query.from_time || Dates.dateHistogramStartDate)
    untilTimeFieldRef.current.value = Dates.FormDate(query.until_time || new Date())
  })

  useEffect(() => {
    const {statistics, metric} = response

    let data = []
    if (!statistics.date_histogram) {
      return
    } else {
      data = Object.keys(statistics.date_histogram).map(key => ({
        time: Dates.JSDate(parseInt(key)),
        value: statistics.date_histogram[key][metric]
      })).filter(d => d.value)
    }

    const fromTime = Dates.JSDate(response.from_time || Dates.dateHistogramStartDate)
    const untilTime = Dates.JSDate(response.until_time || new Date())
    const interval = response.dateHistogramInterval
    const clickable = (interval * Dates.buckets) > 3600

    const handleItemClicked = item => {
      if (!clickable) {
        return
      }
      const fromTime = item.time
      const untilTime = Dates.addSeconds(fromTime, interval)
      setQuery({
        from_time: Dates.APIDate(fromTime),
        until_time: Dates.APIDate(untilTime)
      })
    }

    const width = containerRef.current.offsetWidth
    const height = 250
    const marginRight = 32
    const marginTop = 16
    const marginBottom = 16

    const y = scalePow().range([height - marginBottom, marginTop]).exponent(scale)
    const max = d3.max(data, d => d.value) || 0
    y.domain([0, max])

    const x = scaleTime()
      .domain([Dates.addSeconds(fromTime, -interval), Dates.addSeconds(untilTime, interval)])
      .rangeRound([marginRight, width])

    const container = d3.select(containerRef.current)
    const tooltip = container.select('.' + classes.tooltip)
      .style('opacity', 0)
    const tooltipContent = container.select('.' + classes.tooltipContent)
    const svg = container.select('svg')
      .attr('width', width)
      .attr('height', height)

    const xAxis = d3.axisBottom(x)
    svg.select('.xaxis').remove()
    svg.append('g')
      .attr('transform', `translate(0,${height - marginBottom})`)
      .attr('class', 'xaxis')
      .call(xAxis)

    svg.select('.xlabel').remove()
    svg.append('text')
      .attr('class', 'xlabel')
      .attr('x', width)
      .attr('y', height - 4)
      .attr('dy', '.35em')
      .attr('font-size', '12px')
      .style('text-anchor', 'end')

    const yAxis = d3.axisLeft(y).ticks(Math.min(max, 5), '.0s')
    svg.select('.yaxis').remove()
    svg.append('g')
      .attr('transform', `translate(${marginRight}, 0)`)
      .attr('class', 'yaxis')
      .call(yAxis)

    const {label, shortLabel} = domain.searchMetrics[metric]

    let withData = svg
      .selectAll('.bar').remove().exit()
      .data(data)

    let item = withData.enter()
      .append('g')

    item
      .append('rect')
      .attr('x', d => x(d.time) + 1)
      .attr('y', y(max))
      .attr('width', d => x(Dates.addSeconds(d.time, interval)) - x(d.time) - 2)
      .attr('class', 'background')
      .style('opacity', 0)
      .attr('height', y(0) - y(max))

    item
      .append('rect')
      .attr('class', 'bar')
      .attr('x', d => x(d.time) + 1)
      .attr('y', d => y(d.value))
      .attr('width', d => x(Dates.addSeconds(d.time, interval)) - x(d.time) - 2)
      .attr('height', d => y(0) - y(d.value))
      .style('fill', nomadSecondaryColor.light)

    if (clickable) {
      item
        .style('cursor', 'pointer')
        .on('click', handleItemClicked)
    }

    item
      .on('mouseover', function(d) {
        d3.select(this).select('.background')
          .style('opacity', 0.08)
        if (tooltips) {
          tooltip.transition()
            .duration(200)
            .style('opacity', 1)
          tooltip
            .style('left', x(d.time) + 'px')
            .style('bottom', '24px')
          tooltipContent.html(
            `${d.time.toLocaleDateString()}-${Dates.addSeconds(d.time, interval).toLocaleDateString()} with ${d.value.toLocaleString()} ${shortLabel || label}`)
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
  })

  const handleDatePickerChange = useCallback((event, key) => {
    try {
      const date = new Date(event.target.value).getTime()
      if (date < Dates.JSDate(Dates.dateHistogramStartDate).getTime()) {
        return
      }
      if (date > new Date().getTime()) {
        return
      }
      const value = Dates.APIDate(new Date(event.target.value))
      setQuery({[key]: value})
    } catch (error) {
    }
  }, [setQuery])

  return <Card classes={{root: classes.root}}>
    <CardHeader
      classes={{root: classes.header}}
      title={title}
      titleTypographyProps={{variant: 'body1'}}
      action={(
        <Grid container alignItems='flex-end' spacing={2}>
          <Grid item>
            <TextField
              inputRef={fromTimeFieldRef}
              label="from time"
              type="date"
              defaultValue={Dates.FormDate(query.from_time || Dates.dateHistogramStartDate)}
              onChange={event => handleDatePickerChange(event, 'from_time')}
              InputLabelProps={{
                shrink: true
              }}
            />
          </Grid>
          <Grid item>
            <TextField
              inputRef={untilTimeFieldRef}
              label="until time"
              type="date"
              defaultValue={Dates.FormDate(query.until_time || new Date())}
              onChange={event => handleDatePickerChange(event, 'until_time')}
              InputLabelProps={{
                shrink: true
              }}
            />
          </Grid>
          <Grid item>
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
          </Grid>
        </Grid>
      )}
    />
    <CardContent classes={{root: classes.content}}>
      <div ref={containerRef}>
        <div className={classes.tooltip}>
          <div className={classes.tooltipContent}></div>
        </div>
        <svg />
      </div>
    </CardContent>
  </Card>
}
UploadsHistogram.propTypes = {
  /**
   * An optional title for the chart. If no title is given, the quantity is used.
   */
  title: PropTypes.string,
  /**
   * An optional scale power that is used as the initial scale before the user
   * changes it. Default is 1 (linear scale).
   */
  initialScale: PropTypes.number,
  /**
   * Set to true to enable tooltips for each value.
   */
  tooltips: PropTypes.bool
}
