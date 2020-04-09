import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Select, MenuItem, Card, CardContent, CardHeader } from '@material-ui/core'
import * as d3 from 'd3'
import { scaleBand, scalePow } from 'd3-scale'
import { formatQuantity, nomadPrimaryColor, nomadSecondaryColor } from '../../config.js'
import SearchContext from '../search/SearchContext'

const unprocessedLabel = 'not processed'
const unavailableLabel = 'unavailable'

const _mapping = {
  'energy_total': 'Total energy',
  'energy_total_T0': 'Total energy (0K)',
  'energy_free': 'Free energy',
  'energy_electrostatic': 'Electrostatic',
  'energy_X': 'Exchange',
  'energy_XC': 'Exchange-correlation',
  'energy_sum_eigenvalues': 'Band energy',
  'dos_values': 'DOS',
  'eigenvalues_values': 'Eigenvalues',
  'volumetric_data_values': 'Volumetric data',
  'electronic_kinetic_energy': 'Kinetic energy',
  'total_charge': 'Charge',
  'atom_forces_free': 'Free atomic forces',
  'atom_forces_raw': 'Raw atomic forces',
  'atom_forces_T0': 'Atomic forces (0K)',
  'atom_forces': 'Atomic forces',
  'stress_tensor': 'Stress tensor',
  'thermodynamical_property_heat_capacity_C_v': 'Heat capacity',
  'vibrational_free_energy_at_constant_volume': 'Free energy (const=V)',
  'band_energies': 'Band energies',
  'spin_S2': 'Spin momentum operator',
  'excitation_energies': 'Excitation energies',
  'oscillator_strengths': 'Oscillator strengths',
  'transition_dipole_moments': 'Transition dipole moments'}

function mapKey(key, short = true) {
  let name = key
  const maxLength = 17
  if (key in _mapping) {
    name = _mapping[key]
  }
  if (name.length > maxLength && short) {
    return name.substring(0, maxLength) + '...'
  } else {
    return name
  }
}

class QuantityHistogramUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    title: PropTypes.string.isRequired,
    width: PropTypes.number.isRequired,
    data: PropTypes.object,
    metric: PropTypes.string.isRequired,
    value: PropTypes.string,
    onChanged: PropTypes.func.isRequired,
    defaultScale: PropTypes.number
  }

  static styles = theme => ({
    root: {},
    content: {
      paddingTop: 0
    },
    tooltip: {
      textAlign: 'center',
      position: 'absolute'
    },
    tooltipContent: {
      display: 'inline',
      color: '#fff',
      padding: '4px 8px',
      fontSize: '0.625rem',
      fontFamily: '"Roboto", "Helvetica", "Arial", sans-serif',
      lineHeight: '1.4em',
      borderRadius: '4px',
      backgroundColor: '#616161'
    }
  })

  constructor(props) {
    super(props)
    this.container = React.createRef()
    this.svgEl = React.createRef()
  }

  state = {
    scalePower: this.props.defaultScale || 0.25
  }

  componentDidMount() {
    // TODO this just a workaround for bad layout on initial rendering
    this.updateChart()
    this.updateChart()
  }

  componentDidUpdate() {
    this.updateChart()
  }

  handleItemClicked(item) {
    if (this.props.value === item.key) {
      this.props.onChanged(null)
    } else {
      this.props.onChanged(item.key)
    }
  }

  updateChart() {
    if (!this.props.data) {
      return
    }

    const {classes} = this.props

    const {scalePower} = this.state
    const selected = this.props.value

    const width = this.container.current.offsetWidth
    const left = this.container.current.offsetLeft
    const top = this.container.current.offsetTop
    const height = Object.keys(this.props.data).length * 32

    const data = Object.keys(this.props.data)
      .map(key => ({
        key: key,
        name: mapKey(key),
        value: this.props.data[key][this.props.metric]
      }))

    // keep the data sorting, but put unavailable and not processed to the end
    const unavailableIndex = data.findIndex(d => d.name === unavailableLabel)
    const unprocessedIndex = data.findIndex(d => d.name === unprocessedLabel)
    if (unavailableIndex !== -1) {
      data.push(data.splice(unavailableIndex, 1)[0])
    }
    if (unprocessedIndex !== -1) {
      data.push(data.splice(unprocessedIndex, 1)[0])
    }

    const y = scaleBand().rangeRound([0, height]).padding(0.1)
    const x = scalePow().range([0, width]).exponent(scalePower)

    // we use at least the domain 0..1, because an empty domain causes a weird layout
    const max = d3.max(data, d => d.value) || 1
    x.domain([0, max])
    y.domain(data.map(d => d.name))

    const tooltip = d3.select(this.container.current)
      .append('div')
      .attr('class', classes.tooltip)
      .style('width', width + 'px')
      .style('opacity', 0)

    const tooltipContent = tooltip
      .append('div')
      .attr('class', classes.tooltipContent)

    let svg = d3.select(this.svgEl.current)
    svg.attr('width', width)
    svg.attr('height', height)

    let withData = svg
      .selectAll('.item')
      .data(data, data => data.name)

    withData.exit().remove()

    const rectColor = d => selected === d.key ? nomadPrimaryColor.main : nomadSecondaryColor.light
    const textColor = d => selected === d.key ? '#FFF' : '#000'

    let item = withData.enter()
      .append('g')
      .attr('class', 'item')

    item
      .append('rect')
      .attr('class', 'bar')
      .attr('x', x(0))
      .attr('y', d => y(d.name))
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
      .attr('y', d => y(d.name) + 4)
      .attr('text-anchor', 'start')
      .style('fill', textColor)
      .text(d => d.name)

    item
      .append('text')
      .attr('class', 'value')
      .attr('dy', y.bandwidth())
      .attr('y', d => y(d.name) - 4)
      .attr('x', d => width - 4)
      .attr('text-anchor', 'end')
      .style('fill', textColor)
      .text(d => formatQuantity(d.value))

    item
      .style('cursor', 'pointer')
      .on('click', d => this.handleItemClicked(d))

    item
      .on('mouseover', function(d) {
        tooltip.transition()
          .duration(200)
          .style('opacity', 0.9)
        tooltip
          .style('left', left + 'px')
          .style('top', (y(d.name) + top + 32) + 'px')
        tooltipContent.html(mapKey(d.key, false))
      })
      .on('mouseout', function(d) {
        tooltip.transition()
          .duration(500)
          .style('opacity', 0)
      })

    const t = d3.transition().duration(500)

    item = withData.transition(t)

    item
      .select('.bar')
      .attr('y', d => y(d.name))
      .attr('width', d => x(d.value) - x(0))
      .attr('height', y.bandwidth())
      .style('fill', rectColor)

    item
      .select('.name')
      .text(d => d.name)
      .attr('y', d => y(d.name) + 4)
      .style('fill', textColor)

    item
      .select('.value')
      .text(d => formatQuantity(d.value))
      .attr('y', d => y(d.name) - 4)
      .attr('x', width - 4)
      .style('fill', textColor)
  }

  render() {
    const { classes, title } = this.props

    return (
      <Card classes={{root: classes.root}}>
        <CardHeader
          title={title}
          titleTypographyProps={{variant: 'body1'}}
          action={(
            <Select
              value={this.state.scalePower}
              onChange={(event) => this.setState({scalePower: event.target.value})}
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
          <div ref={this.container}>
            <svg ref={this.svgEl} />
          </div>
        </CardContent>
      </Card>
    )
  }
}

export const QuantityHistogram = withStyles(QuantityHistogramUnstyled.styles)(QuantityHistogramUnstyled)

class QuantityUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    quantity: PropTypes.string.isRequired,
    metric: PropTypes.string.isRequired,
    title: PropTypes.string,
    scale: PropTypes.number
  }
  static styles = theme => ({
    root: {
      marginTop: theme.spacing.unit * 2
    }
  })

  static contextType = SearchContext.type

  render() {
    const {classes, scale, quantity, title, ...props} = this.props
    const {state: {response, query}, setQuery} = this.context

    return <QuantityHistogram
      classes={{root: classes.root}}
      width={300}
      defaultScale={scale || 1}
      title={title || quantity}
      data={response.statistics[quantity]}
      value={query[quantity]}
      onChanged={selection => setQuery({...query, [quantity]: selection})}
      {...props} />
  }
}

export const Quantity = withStyles(QuantityUnstyled.styles)(QuantityUnstyled)
