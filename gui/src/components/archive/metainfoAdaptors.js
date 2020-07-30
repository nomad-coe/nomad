import React, { useMemo, useEffect, useRef } from 'react'
import PropTypes from 'prop-types'
import { useRecoilValue } from 'recoil'
import Adaptor from './adaptors'
import { Item, Content, Compartment, viewConfigState, filterConfigState } from './ArchiveBrowser'
import { Typography, Box, makeStyles } from '@material-ui/core'
import Markdown from '../Markdown'
import { metainfoDef, resolveRef, vicinityGraph } from './metainfo'
import * as d3 from 'd3'

import blue from '@material-ui/core/colors/blue';
import teal from '@material-ui/core/colors/teal';
import lime from '@material-ui/core/colors/lime';
import purple from '@material-ui/core/colors/purple';

export function metainfoAdaptorFactory(obj) {
  if (obj.m_def === 'Section') {
    return new SectionDefAdaptor(obj)
  } else if (obj.m_def === 'SubSection') {
    return new SubSectionDefAdaptor(obj)
  } else if (obj.m_def === 'Quantity') {
    return new QuantityDefAdaptor(obj)
  } else {
    throw new Error('Unknown metainfo definition type')
  }
}

class MetainfoAdaptor extends Adaptor {
  itemAdaptor(key) {
    if (key === '_metainfo') {
      return metainfoAdaptorFactory(metainfoDef(this.e.m_def))
    } else if (this.e[key]) {
      return metainfoAdaptorFactory(resolveRef(this.e[key]))
    } else {
      throw new Error('Unknown item key')
    }
  }
}

export class SectionDefAdaptor extends MetainfoAdaptor {
  itemAdaptor(key) {
    const property = this.e._properties[key]
    if (!property) {
      return super.itemAdaptor(key)
    } else {
      return metainfoAdaptorFactory(property)
    }
  }
  render() {
    return <SectionDef def={this.e} />
  }
}

class QuantityDefAdaptor extends MetainfoAdaptor {
  render() {
    return <QuantityDef def={this.e} />
  }
}

class SubSectionDefAdaptor extends MetainfoAdaptor {
  render() {
    return <SubSectionDef def={this.e} />
  }
}

function SectionDef({def}) {
  const filterConfig = useRecoilValue(filterConfigState)
  const filter = filterConfig.showCodeSpecific ? def => true : def => !def.name.startsWith('x_')
  return <Content style={{backgroundColor: 'grey'}}>
    <Compartment>
      <Definition def={def} isDefinition/>
    </Compartment>
    <Compartment title="visualization">
      <VicinityGraph def={def} />
    </Compartment>
    <Compartment title="sub section definitions">
      {def.sub_sections.filter(filter)
        .map(subSectionDef => {
          const key = subSectionDef.name
          return <Item key={key} itemKey={key}>
            <Typography component="span">
              <Box fontWeight="bold" component="span">
                {subSectionDef.name}
              </Box>
            </Typography>
          </Item>
        })
      }
    </Compartment>
    <Compartment title="quantity definitions">
      {def.quantities.filter(filter)
        .map(quantityDef => {
          const key = quantityDef.name
          return <Item key={key} itemKey={key}>
            <Box component="span" whiteSpace="nowrap">
              <Typography component="span">
                <Box fontWeight="bold" component="span">
                  {quantityDef.name}
                </Box>
              </Typography>
            </Box>
          </Item>
        })
      }
    </Compartment>
  </Content>
}
SectionDef.propTypes = ({
  def: PropTypes.object
})

function SubSectionDef({def}) {
  return <Content>
    <Compartment>
      <Definition def={def} isDefinition/>
    </Compartment>
    <Item itemKey="sub_section">
      <Typography component="span">
        <Box fontWeight="bold" component="span">
          section
        </Box>
      </Typography>: {resolveRef(def.sub_section).name}
    </Item>
  </Content>
}
SubSectionDef.propTypes = ({
  def: PropTypes.object
})

function QuantityDef({def}) {
  return <Content>
    <Compartment>
      <Definition def={def} isDefinition/>
    </Compartment>
    <Compartment title="visualization">
      <VicinityGraph def={def} />
    </Compartment>
  </Content>
}
QuantityDef.propTypes = ({
  def: PropTypes.object
})

const definitionLabels = {
  'Section': 'section',
  'Quantity': 'quantity',
  'SubSection': 'sub section'
}
export function Definition({def, isDefinition}) {
  const viewConfig = useRecoilValue(viewConfigState)
  const color = isDefinition ? 'primary' : 'initial'
  if (viewConfig.showDefinitions) {
    return <Box>
      <Typography color={color} variant="h6">{def.name}</Typography>
      <Typography color={color} variant="caption">{definitionLabels[def.m_def]}{isDefinition ? ' definition' : ''}</Typography>
      {def.description && def.description !== 'None' && viewConfig.showDescriptions &&
        <Box marginTop={1}>
          <Markdown>{def.description}</Markdown>
        </Box>
      }
      <Box marginTop={1}>
        <Item itemKey="_metainfo">
          <Typography>
            metainfo definition
          </Typography>
        </Item>
      </Box>
    </Box>
  } else {
    return ''
  }
}
Definition.propTypes = ({
  def: PropTypes.object,
  isDefinition: PropTypes.bool
})

const useVicinityGraphStyles = makeStyles(theme => ({
  root: {
    '& .links line': {
      // stroke: '#000'
      strokeWidth: 3
    },
    '& text': {
      fontSize: 10
    }
  }
}))
function VicinityGraph({def}) {
  const linkColors = {
    'Quantity': purple[500]
  }
  const nodeColors = {
    'Quantity': lime[500],
    'Section': blue[500],
    'Category': teal[500]
  }
  const graph = useMemo(() => vicinityGraph(def), [def])
  console.log(graph)

  const classes = useVicinityGraphStyles()
  const svgRef = useRef()
  useEffect(() => {
    console.log('#####')
    const svg = d3.select(svgRef.current)
    const width = svg.attr("width")
    const height = svg.attr("height")

    var color = d3.scaleOrdinal(d3.schemeAccent);

    var simulation = d3.forceSimulation()
      .force("link", d3.forceLink().id(d => d.id).distance(75).strength(1))
      .force("charge", d3.forceManyBody().strength(-400))
      .force("center", d3.forceCenter(width / 2, height / 2))
      .force("x", d3.forceX())
      .force("y", d3.forceY())

    var link = svg.select(".links")
      .selectAll("line")
      .data(graph.links)
      .enter().append("line")
      .attr('marker-end','url(#arrowhead)')
      .attr("stroke", d => linkColors[d.def.m_def] || '#000')

    var node = svg.select(".nodes")
      .selectAll("g")
      .data(graph.nodes)
      .enter().append("g")

    node.append("circle")
      .attr("r", 10)
      .attr("fill", d => nodeColors[d.def.m_def] || '#000')
      .call(d3.drag()
        .on("start", dragstarted)
        .on("drag", dragged)
        .on("end", dragended))

    node.append("text")
      .text(d => d.def.name)
      .attr('x', 6)
      .attr('y', 3)

    simulation
      .nodes(graph.nodes)
      .on("tick", ticked)

    simulation.force("link")
      .links(graph.links)

    function ticked() {
      link
        .attr("x1", d => d.source.x)
        .attr("y1", d => d.source.y)
        .attr("x2", d => d.target.x)
        .attr("y2", d => d.target.y)

      node
        .attr("transform", d => "translate(" + d.x + "," + d.y + ")")
    }

    function dragstarted(d) {
      d3.event.sourceEvent.stopPropagation()
      if (!d3.event.active) {
        simulation.alphaTarget(0.3).restart()
      }
      d.fx = d.x
      d.fy = d.y
    }

    function dragged(d) {
      d.fx = d3.event.x
      d.fy = d3.event.y
    }

    function dragended(d) {
      if (!d3.event.active) {
        simulation.alphaTarget(0)
      }
      d.fx = null
      d.fy = null
    }
  }, [graph])
  return <svg className={classes.root} width="300" height="300" ref={svgRef}>
    <g className="links" />
    <g className="nodes" />
    <marker
      id="arrowhead"
      viewBox="-0 -5 10 10"
      refX="21"
      refY="0"
      orient="auto"
      markerWidth="3"
      markerHeight="4"
      xoverflow="visible">
        <path d="M 0,-5 L 10 ,0 L 0,5" fill="#000" stroke="#000" />
      </marker>
  </svg>
}
Definition.propTypes = ({
  def: PropTypes.object.isRequired
})
