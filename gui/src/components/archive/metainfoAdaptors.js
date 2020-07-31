import React, { useMemo, useEffect, useRef, useLayoutEffect } from 'react'
import PropTypes from 'prop-types'
import { useRecoilValue } from 'recoil'
import Adaptor from './adaptors'
import { Item, Content, Compartment, configState } from './ArchiveBrowser'
import { Typography, Box, makeStyles } from '@material-ui/core'
import { metainfoDef, resolveRef, vicinityGraph } from './metainfo'
import * as d3 from 'd3'

import blue from '@material-ui/core/colors/blue';
import teal from '@material-ui/core/colors/teal';
import lime from '@material-ui/core/colors/lime';
import purple from '@material-ui/core/colors/purple';
import grey from '@material-ui/core/colors/grey';
import Markdown from '../Markdown'

export function metainfoAdaptorFactory(obj) {
  if (obj.m_def === 'Section') {
    return new SectionDefAdaptor(obj)
  } else if (obj.m_def === 'SubSection') {
    return new Error('SubSections are not represented in the browser')
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
      if (property.m_def === 'SubSection') {
        return metainfoAdaptorFactory(resolveRef(property.sub_section))
      }
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

function SectionDef({def}) {
  const config = useRecoilValue(configState)
  const filter = config.showCodeSpecific ? def => true : def => !def.name.startsWith('x_')
  return <Content style={{backgroundColor: 'grey'}}>
    <Title def={def} isDefinition/>
    {def.description &&
      <Box marginTop={1} marginBottom={1}>
        <Markdown>{def.description}</Markdown>
      </Box>
    }
    <Compartment>
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
              </Box>{subSectionDef.repeats && <span>&nbsp;(repeats)</span>}
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
    {/* <Meta def={metainfoDef(def.m_def)} /> */}
  </Content>
}
SectionDef.propTypes = ({
  def: PropTypes.object
})

function QuantityDef({def}) {
  return <Content>
    <Title def={def} isDefinition />
    <Compartment title="visualization">
      <VicinityGraph def={def} />
    </Compartment>
    {/* <Meta def={metainfoDef(def.m_def)} /> */}
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
export function Title({def, isDefinition}) {
  const color = isDefinition ? 'primary' : 'initial'
  return <Compartment>
    <Typography color={color} variant="h6">{def.name}</Typography>
    <DefinitionLabel def={def} isDefinition={isDefinition} variant="caption" color={color} />
  </Compartment>
}
Title.propTypes = ({
  def: PropTypes.object.isRequired,
  isDefinition: PropTypes.bool
})

export function DefinitionLabel({def, isDefinition, ...props}) {
  const color = isDefinition ? 'primary' : 'initial'
  return <Typography {...props}>{definitionLabels[def.m_def]}{isDefinition ? ' definition' : ''}</Typography>
}
DefinitionLabel.propTypes = ({
  def: PropTypes.object.isRequired,
  isDefinition: PropTypes.bool
})

const useMetaStyles = makeStyles(theme => ({
  description: {
    marginTop: theme.spacing(1)
  },
  graph: {
    marginTop: theme.spacing(3)
  },
  metainfo: {
    marginBottom: theme.spacing(2)
  },
  metainfoItem: {
    fontWeight: 'bold'
  }
}))
export function Meta({def}) {
  const classes = useMetaStyles()
  const config = useRecoilValue(configState)
  if (!config.showMeta) {
    return ''
  }
  return <Compartment title="meta" color="primary">
    <div className={classes.metainfo}>
      <Item itemKey="_metainfo">
        <DefinitionLabel classes={{root: classes.metainfoItem}} def={def} isDefinition component="span" />
      </Item>
    </div>
    <Markdown classes={{root: classes.description}}>{def.description}</Markdown>
    <div className={classes.graph}>
      <VicinityGraph def={def} />
    </div>
  </Compartment>
}
Meta.propTypes = ({
  def: PropTypes.object
})


const useVicinityGraphStyles = makeStyles(theme => ({
  root: {
    with: '100%',
    minWidth: 300,
    '& .links line': {
      // stroke: '#000'
      strokeWidth: 3
    },
    '& text': {
      fontFamily: 'Titillium Web, sans',
      fontSize: 12,
      fontWeight: 'bold'
    }
  }
}))
export function VicinityGraph({def}) {
  const linkColors = {
    'Quantity': purple[500],
    '_none': grey[700]
  }
  const nodeColors = {
    'Quantity': lime[500],
    'Section': blue[500],
    'Category': teal[500]
  }
  const graph = useMemo(() => {
    const graph = vicinityGraph(def)
    let x1 = Math.min(...graph.nodes.map(n => n.x)) - 32
    let y1 = Math.min(...graph.nodes.map(n => n.y)) - 24
    let x2 = Math.max(...graph.nodes.map(n => n.x)) + 32
    let y2 = Math.max(...graph.nodes.map(n => n.y)) + 24
    const w = Math.max(200, Math.abs(x2 - x1))
    const h = Math.abs(y2 - y1)
    const px = w < 400 ? (400 - w)/2 : 0
    x1 -= px; x2 += px
    graph.viewBox = `${x1} ${y1} ${x2-x1} ${y2-y1}`
    graph.aspectRatio = Math.abs((x2-x1) / (y2-y1))

    console.log(w, h, px, x1, x2)
    return graph
  }, [def])
  console.log(graph)

  const classes = useVicinityGraphStyles()
  const svgRef = useRef()
  useEffect(() => {
    console.log(graph)
    const svg = d3.select(svgRef.current)

    const link = svg.select(".links")
      .selectAll("line")
      .data(graph.links)
      .enter().append("line")
      .attr('marker-end', d => `url(#arrowhead-${d.def.m_def || '_none'})`)
      .attr("stroke", d => linkColors[d.def.m_def || '_none'])

    const node = svg.select(".nodes")
      .selectAll("g")
      .data(graph.nodes)
      .enter().append("g")

    node.append("circle")
      .attr("r", 10)
      .attr("fill", d => nodeColors[d.def.m_def] || '#000')
      .call(d3.drag()
        .on("drag", d => {
          d.x = d3.event.x
          d.y = d3.event.y
          ticked()
        }))

    node.append("text")
      .text(d => d.def.name)
      .attr('text-anchor', 'middle')
      .attr('y', d => ((d.i % 2) === 0) ? 20 : -14)

    function ticked() {
      link
        .attr("x1", d => d.source.x)
        .attr("y1", d => d.source.y)
        .attr("x2", d => d.target.x)
        .attr("y2", d => d.target.y)

      node
        .attr("transform", d => "translate(" + d.x + "," + d.y + ")")
    }
    ticked()
  }, [graph, linkColors, nodeColors])

  useLayoutEffect(() => {
    svgRef.current.style.height = svgRef.current.clientWidth / graph.aspectRatio
  }, [graph, svgRef])

  return <svg
      className={classes.root} ref={svgRef}
      viewBox={graph.viewBox}
  >
    <g className="links" />
    <g className="nodes" />
    {Object.keys(linkColors).map(colorKey => {
      const color = linkColors[colorKey]
      return <marker
        key={colorKey}
        id={`arrowhead-${colorKey}`}
        viewBox="-0 -5 10 10"
        refX="21"
        refY="0"
        orient="auto"
        markerWidth="3"
        markerHeight="4"
        xoverflow="visible">
          <path d="M 0,-5 L 10 ,0 L 0,5" fill={color} stroke={color} />
      </marker>
    })}
  </svg>
}
VicinityGraph.propTypes = ({
  def: PropTypes.object.isRequired
})
