import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Select, MenuItem } from '@material-ui/core'
import Grid from '@material-ui/core/Grid'
import TextField from '@material-ui/core/TextField'
import * as d3 from 'd3'
import { scaleBand, scalePow } from 'd3-scale'
import { nomadSecondaryColor } from '../../config.js'
import SearchContext from './SearchContext'
import { compose } from 'recompose'
import { withApi } from '../api'
import { Quantity } from './QuantityHistogram'

class UploadsHistogramUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    height: PropTypes.number.isRequired,
    data: PropTypes.object,
    interval: PropTypes.string,
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
      interval: this.props.interval || '1M',
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

  intervals = [
    {
      label: 'Yearly',
      value: '1y',
      number: 31536000000
    },
    {
      label: 'Monthly',
      value: '1M',
      number: 2678400000
    },
    {
      label: 'Daily',
      value: '1d',
      number: 86400000
    },
    {
      label: 'Hourly',
      value: '1h',
      number: 3600000
    },
    {
      label: 'Minute',
      value: '1m',
      number: 60000
    },
    {
      label: 'Second',
      value: '1s',
      number: 1000
    }
  ]

  timeInterval = Object.assign({}, ...this.intervals.map(e => ({[e.value]: e.number})))

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
    const interval = this.state.interval
    const from_time = new Date(this.state.from_time)
    const until_time = new Date(this.state.until_time)
    this.props.onChanged(from_time.toISOString(), until_time.toISOString(), interval)
  }

  handleIntervalChange(newInterval) {
    // TODO: add a refresh button so directly updating interval is not necessary
    this.setState({interval: newInterval}, () => this.handleQueryChange())
  }

  handleTimeChange(newTime, key, target) {
    let date
    if (!newTime) {
      date = key === 'from_time' ? new Date(this.startDate) : new Date()
    } else {
      date = new Date(newTime)
    }
    if (target === 'state' || target === 'all') {
      key === 'from_time' ? this.setState({from_time: date.getTime()}) : this.setState({until_time: date.getTime()})
    }
    if (target === 'picker' || target === 'all') {
      document.getElementById(key).value = date.toISOString().substring(0, 10)
    }
  }

  handleItemClicked(item) {
    const selected = item.time
    if (selected === this.state.time) {
      this.props.onChanged(null, null, null)
    } else {
      const deltaT = this.timeInterval[this.state.interval]
      this.handleTimeChange(selected, 'from_time', 'all')
      this.handleTimeChange(selected + deltaT, 'until_time', 'all')
      this.handleQueryChange()
    }
  }

  resolveDate(name) {
    const date = new Date(parseInt(name, 10))
    const year = date.toLocaleDateString(undefined, {year: 'numeric'})
    const month = date.toLocaleDateString(undefined, {month: 'short'})
    const day = date.toLocaleDateString(undefined, {day: 'numeric'})
    const hour = date.toLocaleTimeString(undefined, {hour: 'numeric'})
    const min = date.toLocaleTimeString(undefined, {minute: 'numeric'})
    const sec = date.toLocaleTimeString(undefined, {second: 'numeric'})

    const intervals = {
      '1y': year,
      '1M': month,
      '1d': day,
      '1h': hour,
      '1m': min,
      '1s': sec
    }

    return intervals[this.state.interval]
  }

  hover(svg, bar) {
    const textOffset = 25

    const tooltip = svg.append('g')
      .attr('class', 'tooltip')
      .style('display', 'none')

    const hoverBox = tooltip.append('rect')
      .attr('width', 10)
      .attr('height', 20)
      .attr('fill', 'white')
      .style('opacity', 0.0)

    const text = tooltip.append('text')
      .attr('x', textOffset)
      .attr('dy', '1.2em')
      .style('text-anchor', 'start')
      .attr('font-size', '12px')
    // .attr('font-weight', 'bold')

    bar
      .on('mouseover', () => {
        tooltip.style('display', null)
        let { width } = text.node().getBBox()
        hoverBox.attr('width', `${width + textOffset}px`)
      })
      .on('mouseout', () => tooltip.style('display', 'none'))
      .on('mousemove', function(d) {
        let xPosition = d3.mouse(this)[0] - 15
        let yPosition = d3.mouse(this)[1] - 25

        tooltip.attr('transform', `translate( ${xPosition}, ${yPosition})`)
        tooltip.attr('data-html', 'true')
        tooltip.select('text').text(new Date(d.time).toISOString() + ': ' + d.value)
      })
  }

  updateChart() {
    let data = []
    if (!this.props.data) {
      return
    } else {
      data = Object.keys(this.props.data).map(key => ({
        time: parseInt(key, 10),
        name: this.resolveDate(key),
        value: this.props.data[key][this.props.metric]
      }))
    }

    data.sort((a, b) => d3.ascending(a.time, b.time))
    if (data.length > 0) {
      this.handleTimeChange(this.state.from_time, 'from_time', 'picker')
      this.handleTimeChange(this.state.until_time, 'until_time', 'picker')
    }

    const scalePower = this.state.scalePower
    const width = this.container.current.offsetWidth
    const height = this.props.height
    const margin = Math.round(0.1 * height)

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

    const yAxis = d3.axisLeft(y)
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
      .attr('y', 0)
      .attr('dy', '1em')
      .attr('text-anchor', 'start')
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
      .on('click', d => this.handleItemClicked(d))

    svg.select('.tooltip').remove()
    svg.call(this.hover, item)
  }

  render() {
    return (
      <div>
        <Grid container justify='space-around' spacing={24}>
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
          <Grid item xs={3}>
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
          <Grid item xs={3}>
            <Select
              id='interval'
              value={this.state.interval}
              onChange={(event) => this.handleIntervalChange(event.target.value)}
              label= 'interval'
            > {this.intervals.map(item => (
                <MenuItem value={item.value} key={item.value}>
                  {item.label}
                </MenuItem>))}
            </Select>
          </Grid>
          <Grid item xs={3}>
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
        </Grid>
        <div ref={this.container}>
          <svg ref={this.svgEl}></svg>
        </div>
      </div>
    )
  }
}

export const UploadsHistogram = withStyles(UploadsHistogramUnstyled.styles)(UploadsHistogramUnstyled)

class UploadersListUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired
  }

  static styles = theme => ({
    root: {
      marginTop: theme.spacing.unit * 2
    }
  })

  static contextType = SearchContext.type

  render() {
    const {state: {usedMetric}} = this.context

    return (
      <Grid>
        <Quantity quantity="uploader" title="Top Uploaders" scale={1} metric={usedMetric} />
      </Grid>
    )
  }
}

export const UploadersList = withStyles(UploadersListUnstyled.styles)(UploadersListUnstyled)

class UploadsChart extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    metricsDefinitions: PropTypes.object.isRequired
  }
  static styles = theme => ({
    root: {
      marginTop: theme.spacing.unit
    }
  })

  static contextType = SearchContext.type

  render() {
    const {classes, metricsDefinitions, ...props} = this.props
    const {state: {response, usedMetric, query}, setQuery} = this.context

    return (
      <Grid container spacing={24}>
        <Grid item xs={12}>
          <UploadsHistogram
            classes={{root: classes.root}}
            height={250}
            defaultScale={1}
            data={response.statistics.date_histogram}
            metric={usedMetric}
            metricsDefinitions={metricsDefinitions}
            interval={'1M'}
            onChanged={(from_time, until_time, interval) => setQuery({...query, from_time: from_time, until_time: until_time, interval: interval})}
            {...props} />
        </Grid>
        <Grid item xs={12}>
          <UploadersList />
        </Grid>
      </Grid>
    )
  }
}

export default compose(withApi(false, false), withStyles(UploadsChart.styles))(UploadsChart)
