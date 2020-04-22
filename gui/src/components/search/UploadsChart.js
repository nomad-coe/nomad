import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Select, MenuItem } from '@material-ui/core'
import Button from '@material-ui/core/Button'
import RefreshIcon from '@material-ui/icons/Refresh'
import Grid from '@material-ui/core/Grid'
import TextField from '@material-ui/core/TextField'
import * as d3 from 'd3'
import { scaleBand, scalePow } from 'd3-scale'
import { nomadSecondaryColor } from '../../config.js'
import { searchContext } from './SearchContext'
import { compose } from 'recompose'
import { withApi } from '../api'

class UploadsHistogramUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    height: PropTypes.number.isRequired,
    data: PropTypes.object,
    metric: PropTypes.string.isRequired,
    metricsDefinitions: PropTypes.object.isRequired,
    onChanged: PropTypes.func.isRequired,
    defaultScale: PropTypes.number
  }

  static styles = theme => ({
    root: {},
    content: {
      paddingTop: 10
    }
  })

  constructor(props) {
    super(props)
    this.state = {
      scalePower: this.props.defaultScale || 1.0,
      time: null,
      from_time: 0,
      until_time: 0
    }

    this.container = React.createRef()
    this.svgEl = React.createRef()
  }

  startDate = '2013-01-01'

  scales = [
    {
      label: 'Linear',
      value: 1.0
    },
    {
      label: '1/2',
      value: 0.5
    },
    {
      label: '1/4',
      value: 0.25
    },
    {
      label: '1/8',
      value: 0.25
    }
  ]

  componentDidMount() {
    const from_time = new Date(this.startDate).getTime()
    const until_time = new Date().getTime()
    this.handleTimeChange(from_time, 'from_time', 'all')
    this.handleTimeChange(until_time, 'until_time', 'all')
  }

  componentDidUpdate() {
    this.updateChart()
  }

  handleQueryChange() {
    const from_time = new Date(this.state.from_time)
    const until_time = new Date(this.state.until_time)
    this.props.onChanged(from_time.toISOString(), until_time.toISOString())
  }

  handleTimeChange(newTime, key, target) {
    let date
    if (!newTime) {
      date = key === 'from_time' ? new Date(this.startDate) : new Date()
    } else {
      date = new Date(newTime)
    }

    if (target === 'state' || target === 'all') {
      if (key === 'from_time') {
        if (this.state.from_time !== date.getTime()) {
          this.setState({from_time: date.getTime()})
        }
      } else if (key === 'until_time') {
        if (this.state.until_time !== date.getTime()) {
          this.setState({until_time: date.getTime()})
        }
      }
    }

    if (target === 'picker' || target === 'all') {
      document.getElementById(key).value = date.toISOString().substring(0, 10)
    }
  }

  handleItemClicked(item, deltaT) {
    const selected = item.time
    if (selected === this.state.time) {
      this.props.onChanged(null, null)
    } else {
      this.handleTimeChange(selected, 'from_time', 'all')
      this.handleTimeChange(selected + deltaT, 'until_time', 'all')
      this.handleQueryChange()
    }
  }

  resolveDate(name, deltaT) {
    const date = new Date(parseInt(name, 10))
    const year = date.toLocaleDateString(undefined, {year: 'numeric'})
    const quarter = Math.floor((date.getMonth() + 3) / 3)
    const month = date.toLocaleDateString(undefined, {month: 'short'})
    const week = (date) => {
      const first = new Date(date.getFullYear(), 0, 1)
      return Math.ceil((((date - first) / 86400000) + first.getDay() + 1) / 7)
    }
    const day = date.toLocaleDateString(undefined, {day: 'numeric'})
    const hour = date.toLocaleTimeString(undefined, {hour: 'numeric'})
    const min = date.toLocaleTimeString(undefined, {minute: 'numeric'})
    const sec = date.toLocaleTimeString(undefined, {second: 'numeric'})

    const times = [31536000, 7776000, 2419200, 604800, 864000, 3600, 60, 1]
    const diffs = times.map(t => Math.abs(t - (deltaT / 1000)))

    const intervals = [year, 'Q' + quarter, month, 'W' + week, day, hour, min, sec]
    return intervals[diffs.indexOf(Math.min(...diffs))]
  }

  updateChart() {
    let data = []
    if (!this.props.data) {
      return
    } else {
      data = Object.keys(this.props.data).map(key => ({
        time: parseInt(key, 10),
        // name: this.resolveDate(key),
        value: this.props.data[key][this.props.metric]
      }))
    }

    data.sort((a, b) => d3.ascending(a.time, b.time))
    if (data.length > 0) {
      this.handleTimeChange(this.state.from_time, 'from_time', 'picker')
      this.handleTimeChange(this.state.until_time, 'until_time', 'picker')
    }

    let deltaT = 31536000000
    if (data.length > 1) {
      deltaT = data[1].time - data[0].time
    }

    data.forEach(d => { d.name = this.resolveDate(d.time, deltaT) })

    const scalePower = this.state.scalePower
    const width = this.container.current.offsetWidth
    const height = this.props.height
    const margin = Math.round(0.15 * height)

    const x = scaleBand().rangeRound([margin, width]).padding(0.1)
    const y = scalePow().range([height - margin, margin]).exponent(scalePower)

    const max = d3.max(data, d => d.value) || 0
    x.domain(data.map(d => d.name))
    y.domain([0, max])

    let svg = d3.select(this.svgEl.current)
    svg.attr('width', width)
    svg.attr('height', height)

    const xAxis = d3.axisBottom(x)
    svg.select('.xaxis').remove()
    svg.append('g')
      .attr('transform', `translate(0,${height - margin})`)
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
      .attr('transform', `translate(${margin}, 0)`)
      .attr('class', 'yaxis')
      .call(yAxis)

    const {label, shortLabel} = this.props.metricsDefinitions[this.props.metric]
    svg.select('.ylabel').remove()
    svg.append('text')
      .attr('class', 'ylabel')
      .attr('x', 0)
      .attr('y', 10)
      .attr('dy', '1em')
      .attr('font-size', '12px')
      .text(`${shortLabel || label}`)

    let withData = svg
      .selectAll('.bar').remove().exit()
      .data(data)

    let item = withData.enter()
      .append('g')

    item
      .append('rect')
      .attr('class', 'bar')
      .attr('x', d => x(d.name))
      .attr('y', d => y(d.value))
      .attr('width', x.bandwidth())
      .attr('height', d => y(0) - y(d.value))
      .style('fill', nomadSecondaryColor.light)

    item
      .style('cursor', 'pointer')
      .on('click', d => this.handleItemClicked(d, deltaT))

    svg.select('.tooltip').remove()
    let tooltip = svg.append('g')
    tooltip
      .attr('class', 'tooltip')
      .style('visibility', 'hidden')

    tooltip.append('rect')
      .attr('x', 0)
      .attr('rx', 6)
      .attr('ry', 6)
      .attr('width', 100)
      .attr('height', 40)
      .attr('fill', 'grey')
      .style('opacity', 1.0)

    let tooltipText = tooltip.append('text')
      .attr('dy', '1.2em')
      .attr('font-family', 'Arial, Helvetica, sans-serif')
      .attr('font-size', '10px')
      .attr('fill', 'white')
      .style('text-anchor', 'middle')

    item
      .on('mouseover', function(d) {
        tooltip.style('visibility', 'visible')
        const date = new Date(d.time)
        const value = `${date.toLocaleDateString()}\n
          ${date.toLocaleTimeString()}
          ${d.value} ${shortLabel || label}`
        tooltipText.selectAll('tspan')
          .data((value).split(/\n/)).join('tspan')
          .attr('x', 50)
          .attr('y', (d, i) => `${i * 1.2}em`)
          .text(d => d)
        const xPosition = x(d.name) + 20
        const yPosition = y(d.value) - 20
        tooltip.attr('transform', `translate( ${xPosition}, ${yPosition})`)
      })
      .on('mouseout', () => tooltip.style('visibility', 'hidden'))
  }

  render() {
    return (
      <div>
        <Grid container justify='space-between' alignItems='flex-end'>
          <Grid item xs={2}>
            <Select
              margin='none'
              id='scales'
              value={this.state.scalePower}
              onChange={(event) => this.setState({scalePower: event.target.value})}
              label= 'scale'
            > {this.scales.map(item => (
                <MenuItem
                  value={item.value}
                  key={item.label}> {item.label}
                </MenuItem>))}
            </Select>
          </Grid>
          <Grid item xs={2}>
            <TextField
              id='from_time'
              label="from time"
              type="date"
              onChange={(event) => this.handleTimeChange(event.target.value, 'from_time', 'state')}
              InputLabelProps={{
                shrink: true
              }}
            />
          </Grid>
          <Grid item xs={2}>
            <TextField
              id='until_time'
              label="until time"
              type="date"
              onChange={(event) => this.handleTimeChange(event.target.value, 'until_time', 'state')}
              InputLabelProps={{
                shrink: true
              }}
            />
          </Grid>
          <Grid item xs={2}>
            <Button
              variant='outlined'
              color='default'
              onClick={() => this.handleQueryChange()}
            > Refresh <RefreshIcon/>
            </Button>
          </Grid>
        </Grid>
        <div ref={this.container}>
          <svg ref={this.svgEl}></svg>
        </div>
      </div>
    )
  }
}

export const UploadsHistogram = withStyles(UploadsHistogramUnstyled.styles)(UploadsHistogramUnstyled)

class UploadsChart extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    metricsDefinitions: PropTypes.object.isRequired
  }
  static styles = theme => ({
    root: {
      marginTop: theme.spacing(1)
    }
  })

  static contextType = searchContext

  componentDidMount() {
    const {setStatisticsToRefresh} = this.context
    setStatisticsToRefresh('date_histogram')
  }

  render() {
    const {classes, metricsDefinitions, ...props} = this.props
    const {response: {statistics, metric}, query, setQuery} = this.context

    return (
      <Grid container spacing={2}>
        <Grid item xs={12}>
          <UploadsHistogram
            classes={{root: classes.root}}
            height={250}
            defaultScale={1}
            data={statistics.date_histogram}
            metric={metric}
            metricsDefinitions={metricsDefinitions}
            onChanged={(from_time, until_time) => setQuery({...query, from_time: from_time, until_time: until_time})}
            {...props} />
        </Grid>
      </Grid>
    )
  }
}

export default compose(withApi(false, false), withStyles(UploadsChart.styles))(UploadsChart)
