import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Typography } from '@material-ui/core'
import * as d3 from 'd3'
import { scaleBand, scaleLinear } from 'd3-scale'
import chroma from 'chroma-js'
import repoColor from '@material-ui/core/colors/deepPurple'

class QuantityHistogram extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    title: PropTypes.string.isRequired,
    width: PropTypes.number.isRequired,
    data: PropTypes.object,
    onSelectionChanged: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {

    }
  })

  constructor(props) {
    super(props)
    this.container = React.createRef()
    this.svgEl = React.createRef()
  }

  state = {
    selected: undefined
  }

  componentDidMount() {
    this.updateChart()
  }

  componentDidUpdate() {
    this.updateChart()
  }

  handleItemClicked(item) {
    const isSelected = this.state.selected === item.name
    let selected
    if (isSelected) {
      selected = undefined
    } else {
      selected = item.name
    }

    this.setState({selected: selected})
    this.props.onSelectionChanged(selected)
  }

  updateChart() {
    if (!this.props.data) {
      return
    }

    const { selected } = this.state

    const width = this.container.current.offsetWidth
    const height = Object.keys(this.props.data).length * 32

    const data = Object.keys(this.props.data).map(key => ({
      name: key,
      value: this.props.data[key]
    }))

    const y = scaleBand().rangeRound([0, height]).padding(0.1)
    const x = scaleLinear().range([0, width])
    const heatmapScale = chroma.scale(['#ffcdd2', '#d50000'])

    x.domain([0, d3.max(data, d => d.value)])
    y.domain(data.map(d => d.name))
    heatmapScale.domain([0, d3.max(data, d => d.value)], 10, 'log')

    let svg = d3.select(this.svgEl.current)
    svg.attr('width', width)
    svg.attr('height', height)

    let withData = svg
      .selectAll('g')
      .data(data, data => data.name)

    withData.exit().remove()

    const rectColor = d => selected === d.name ? repoColor[500] : heatmapScale(d.value)
    const textColor = d => selected === d.name ? '#FFF' : '#000'

    let item = withData.enter()
      .append('g')

    item
      .append('rect')
      .attr('x', x(0))
      .attr('y', d => y(d.name))
      .attr('width', d => x(d.value) - x(0))
      .attr('height', y.bandwidth())
      .style('fill', rectColor)
      .style('stroke', '#000')
      .style('stroke-width', '1px')
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
      .attr('x', d => x(d.value) - 4)
      .attr('text-anchor', 'end')
      .style('fill', textColor)
      .text(d => '' + d.value)

    item
      .style('cursor', 'pointer')
      .on('click', d => this.handleItemClicked(d))

    const t = d3.transition().duration(500)

    item = withData.transition(t)

    item
      .select('rect')
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
      .text(d => '' + d.value)
      .attr('y', d => y(d.name) - 4)
      .attr('x', d => x(d.value) - 4)
      .style('fill', textColor)
  }

  render() {
    const { classes, title } = this.props

    return (
      <div className={classes.root} ref={this.container}>
        <Typography variant="body1">{title}</Typography>
        <svg ref={this.svgEl} />
      </div>
    )
  }
}

export default withStyles(QuantityHistogram.styles)(QuantityHistogram)
