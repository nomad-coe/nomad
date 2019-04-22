import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Select, MenuItem, Card, CardContent, CardHeader } from '@material-ui/core'
import * as d3 from 'd3'
import { scaleBand, scalePow } from 'd3-scale'
import chroma from 'chroma-js'
import repoColor from '@material-ui/core/colors/deepPurple'
import { formatQuantity } from '../../config.js'

const unprocessed_label = 'not processed'
const unavailable_label = 'unavailable'

class QuantityHistogram extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    title: PropTypes.string.isRequired,
    width: PropTypes.number.isRequired,
    data: PropTypes.object,
    metric: PropTypes.string.isRequired,
    value: PropTypes.string,
    onChanged: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {},
    content: {
      paddingTop: 0
    }
  })

  constructor(props) {
    super(props)
    this.container = React.createRef()
    this.svgEl = React.createRef()
  }

  state = {
    scalePower: 0.25
  }

  componentDidMount() {
    this.updateChart()
  }

  componentDidUpdate() {
    this.updateChart()
  }

  handleItemClicked(item) {
    if (this.props.value === item.name) {
      this.props.onChanged(null)
    } else {
      this.props.onChanged(item.name)
    }
  }

  updateChart() {
    if (!this.props.data) {
      return
    }

    const { scalePower } = this.state
    const selected = this.props.value

    const width = this.container.current.offsetWidth
    const height = Object.keys(this.props.data).length * 32

    const data = Object.keys(this.props.data).map(key => ({
      name: key,
      value: this.props.data[key][this.props.metric]
    }))

    data.sort((a, b) => {
      const nameA = a.name
      const nameB = b.name

      if (nameA === nameB) {
        return 0
      }

      if (nameA === unprocessed_label) {
        return 1
      }
      if (nameB === unprocessed_label) {
        return -1
      }
      if (nameA === unavailable_label) {
        return 1
      }
      if (nameB === unavailable_label) {
        return -1
      }
      return nameA.localeCompare(nameB)
    })

    const y = scaleBand().rangeRound([0, height]).padding(0.1)
    const x = scalePow().range([0, width]).exponent(scalePower)
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
      .text(d => formatQuantity(d.value))

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

export default withStyles(QuantityHistogram.styles)(QuantityHistogram)