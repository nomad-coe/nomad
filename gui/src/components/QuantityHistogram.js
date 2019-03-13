import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core'
import * as d3 from 'd3'
import { scaleBand, scaleLinear } from 'd3-scale'
import chroma from 'chroma-js'

class QuantityHistogram extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    width: PropTypes.number.isRequired,
    height: PropTypes.number.isRequired,
    data: PropTypes.object.isRequired
  }

  static styles = theme => ({
    root: {

    }
  })

  constructor(props) {
    super(props)
    this.svgEl = React.createRef()
  }

  componentDidMount() {
    this.updateChart()
  }

  componentDidUpdate() {
    this.updateChart()
  }

  updateChart() {
    const { width, height } = this.props

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

    let withData = svg
      .selectAll('g')
      .data(data, data => data.name)

    withData.exit().remove()

    let item = withData.enter()
      .append('g')

    item
      .append('rect')
      .attr('x', x(0))
      .attr('y', d => y(d.name))
      .attr('width', d => x(d.value) - x(0))
      .attr('height', y.bandwidth())
      .style('fill', d => heatmapScale(d.value))

    item
      .append('text')
      .attr('class', 'name')
      .attr('dy', '.75em')
      .attr('x', x(0) + 4)
      .attr('y', d => y(d.name) + 4)
      .attr('text-anchor', 'start')
      .text(d => d.name)

    item
      .append('text')
      .attr('class', 'value')
      .attr('dy', y.bandwidth())
      .attr('y', d => y(d.name) - 4)
      .attr('x', d => x(d.value) - 4)
      .attr('text-anchor', 'end')
      .text(d => '' + d.value)

    const t = d3.transition().duration(500)

    item = withData.transition(t)

    item
      .select('rect')
      .attr('y', d => y(d.name))
      .attr('width', d => x(d.value) - x(0))
      .attr('height', y.bandwidth())
      .style('fill', d => heatmapScale(d.value))

    item
      .select('.name')
      .text(d => d.name)
      .attr('y', d => y(d.name) + 4)

    item
      .select('.value')
      .text(d => '' + d.value)
      .attr('y', d => y(d.name) - 4)
      .attr('x', d => x(d.value) - 4)
  }

  render() {
    const { classes } = this.props

    return (
      <div className={classes.root}>
        <svg
          width={this.props.width}
          height={this.props.height}
          ref={this.svgEl} >
        </svg>
      </div>
    )
  }
}

export default withStyles(QuantityHistogram.styles)(QuantityHistogram)
