/*
 * Copyright The NOMAD Authors.
 *
 * This file is part of NOMAD. See https://nomad-lab.eu for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License")
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
import React, { useEffect, useMemo, useRef, useState } from 'react'
import PropTypes from 'prop-types'
import { makeStyles, Tooltip, IconButton } from '@material-ui/core'
import grey from '@material-ui/core/colors/grey'
import red from '@material-ui/core/colors/red'
import { PropertyCard, PropertyGrid } from './PropertyCard'
import { resolveNomadUrl, resolveInternalRef, createEntryUrl } from '../../../utils'
import { useApi } from '../../api'
import { useHistory } from 'react-router-dom'
import * as d3 from 'd3'
import { getUrl } from '../../nav/Routes'
import { nomadFontFamily, nomadPrimaryColor, nomadSecondaryColor, apiBase } from '../../../config'
import { Replay, Undo, Redo } from '@material-ui/icons'

const useWorkflowGraphStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    minWidth: 1000,
    '& .link': {
      strokeWidth: 3,
      fill: 'none',
      stroke: grey[500],
      strokeOpacity: 0.5
    },
    '& .crosslink': {
      fill: 'none',
      stroke: red[500],
      strokeOpacity: 0.5,
      strokeWidth: 3
    },
    '& .text': {
      fontFamily: nomadFontFamily,
      fontSize: 12,
      fontWeight: 'bold',
      dy: '0.35em',
      fill: grey[800]
    },
    '& .circle': {
      stroke: grey[500],
      strokeWidth: 2,
      strokeOpacity: 0.5
    }
  }
}))

function addMarkers(svg, nodes, offset) {
  const markerWidth = 5

  const markerIds = []

  const addIds = (node) => {
    const nodeLinks = node.data.crossLinks || []
    nodeLinks.forEach(link => {
      const nodeTarget = nodes.filter(node => node.data.key === link[0])
      nodeTarget.forEach(target => {
        markerIds.push(`${node.id}-${target.id}`)
      })
    })
    if (node._children) {
      node._children.forEach(child => addIds(child))
    }
  }

  nodes.forEach(node => {
    addIds(node)
  })

  svg.append('defs').selectAll('marker')
    .data(markerIds)
    .enter()
    .append('marker')
    .attr('class', d => `marker-${d}`)
    .attr('id', String)
    .attr('viewBox', '0 -5 10 10')
    .attr('refX', offset || 0)
    .attr('refY', 0)
    .attr('markerWidth', markerWidth)
    .attr('markerHeight', markerWidth)
    .attr('xoverflow', 'visible')
    .attr('orient', 'auto')
    .style('fill', red[500])
    .style('fill-opacity', 0.5)
    .append('path')
    .attr('d', 'M0,-5L10,0L0,5')
}

const ForceDirected = React.memo(({data, layout, setTooltipContent}) => {
  const classes = useWorkflowGraphStyles()
  const svgRef = useRef()
  const history = useHistory()
  const finalLayout = useMemo(() => {
    const defaultLayout = {
     width: 600,
     margin: {top: 40, bottom: 40, left: 20, right: 20},
     circleRadius: 17,
     linkDistance: 40
    }
    return {...defaultLayout, ...layout}
  }, [layout])

  useEffect(() => {
    const svg = d3.select(svgRef.current)
    svg.selectAll('g').remove()
    const { width, circleRadius, linkDistance, margin } = finalLayout
    svg
      .attr('width', width)
      .attr('height', width)

    const svgGroup = svg.append('g')

    const gCrossLink = svgGroup.append('g')
      .attr('class', 'crosslink')

    const gLink = svgGroup.append('g')
      .attr('class', 'link')

    const gNode = svgGroup.append('g')
      .attr('class', 'node')
      .attr('cursor', 'pointer')
      .attr('pointer-events', 'all')

    // fix inputs and outputs to edge of frame
    const mid = width / 2
    const dy = circleRadius * 10

    if (!data.children) data.children = []

    const inputChildren = data.children.filter(d => d.intent && d.intent.startsWith('input'))
    let offset = inputChildren.length * dy / 2
    let fixPoints = inputChildren.map((child, index) => mid + index * dy - offset)
    for (const [index, child] of inputChildren.entries()) {
      // vertical configuration, switch fixX and fixY for vertical
      child.fixX = fixPoints[index]
      child.fixY = margin.top
    }
    const outputChildren = data.children.filter(d => d.intent && d.intent.startsWith('output'))
    offset = outputChildren.length * dy / 2
    fixPoints = outputChildren.map((child, index) => mid + index * dy - offset)
    for (const [index, child] of outputChildren.entries()) {
      // vertical configuration, switch fixX and fixY for vertical
      child.fixX = fixPoints[index]
      child.fixY = width - margin.bottom
    }

    let nodes = []
    function flatten(node) {
      node.children = node.children || []
      node.children.forEach(child => flatten(child))
      nodes.push(node)
    }

    flatten(data)

    // add crosslinks between nodes
    function addCrossLinksData(data) {
      const workflows = []
      const tasks = []
      const inputs = []
      const outputs = []
      data.children.forEach(child => {
        if (child.type.startsWith('workflow')) workflows.push(child)
        else if (child.type.startsWith('task')) tasks.push(child)
        else if (child.intent && child.intent.startsWith('input')) inputs.push(child)
        else if (child.intent && child.intent.startsWith('output')) outputs.push(child)
      })

      data.children = [...inputs, ...outputs, ...workflows, ...tasks]

      const addLink = (link, node, dash = '3,3') => {
        link.push(dash)
        nodes.forEach(d => {
          d.crossLinks = d.crossLinks || []
          if (d.key === node.key) {
            // apply link to duplicate nodes
            const links = d.crossLinks.map(l => l[0])
            if (link[0] && !links.includes(link[0])) d.crossLinks.push(link)
          }
          if (d.key === link[0]) {
            node.crossLinks = node.crossLinks || []
            const links = node.crossLinks.map(l => l[0])
            if (link[0] && !links.includes(link[0])) node.crossLinks.push(link)
          }
        })
      }

      const allTasks = [...workflows, ...tasks]

      // add link from inputs to workflow
      inputs.forEach(input => addLink([data.key, 'Input'], input))
      // add link from workflow to outputs
      outputs.forEach(output => addLink([output.key, 'Output'], data))
      // add link from from input to first task
      if (tasks.length > 0) {
        // add link from input to first task
        inputs.forEach(input => {
          addLink([tasks[0].key, 'Input'], input)
        })
        // add link from last task to output
        outputs.forEach(output => {
          // add link only if no task is connected to output
          let linked = false
          for (const node of allTasks) {
            const links = node.crossLinks.map(l => l[0])
            if (links.includes(output.key)) {
              linked = true
              break
            }
          }
          if (!linked) addLink([output.key, 'Output'], tasks[tasks.length - 1])
        })
      }

      // add links between tasks if inputs/outputs are connected
      allTasks.forEach(task1 => {
        allTasks.forEach(task2 => {
          if (task1.key !== task2.key) {
            const outputs = task1.children
              .filter(child => child.intent && child.intent.startsWith('output'))
              .map(child => child.key)
            const inputs = task2.children
              .filter(child => child.intent && child.intent.startsWith('input'))
              .map(child => child.key)
            let linked = false
            for (const input of inputs) {
              if (outputs.includes(input)) {
                linked = true
                break
              }
            }
            if (linked) addLink([task2.key, 'Sequential tasks'], task1)
          }
        })
      })

      // iterate over all tasks
      workflows.forEach(workflow => addCrossLinksData(workflow))
      tasks.forEach(task => addCrossLinksData(task))
    }

    addCrossLinksData(data)

    // pack layout for the task quantities
    const pack = d3.pack()
      .size([circleRadius * 2, circleRadius * 2])
      .padding(5)

    // define transition of nodes
    // const transition = svg.transition()
    //   .duration(250)

    function addSubGraph(node) {
      if (!node._children || !node._children.length) return

      const d = {...node.data}
      d.children = (node._children || []).map(child => {
        return {...child.data, crossLinks: []}
      })
      d.crossLinks = []
      d.children.forEach(child => {
        if (child.intent === 'output') d.crossLinks.push([`${child.key}_subgraph`, ''])
        if (child.intent === 'input') child.crossLinks.push([`${d.key}_subgraph`, ''])
      })

      const dRoot = d3.hierarchy(d)
      dRoot.sum(d => d.size).sort((a, b) => a.data.index - b.data.index)
      dRoot.descendants().forEach((node, i) => {
        node.data.key = `${node.data.key}_subgraph`
        node.id = i + nodes.length + subGraphNodes.length
      })
      pack(dRoot)

      const crossLinks = []
      dRoot.children = dRoot.children || []
      dRoot.children.forEach(child => {
        let label = ''
        if (child.data.intent === 'output') {
          for (const link of node.data.crossLinks) {
            if (child.data.key.startsWith(link[0])) {
              label = link[1]
              break
            }
          }
          crossLinks.push({source: dRoot, target: child, id: `${dRoot.id}-${child.id}`, label: label})
        }
        if (child.data.intent === 'input') {
          let label = ''
          for (const nodeChild of node._children) {
            for (const link of nodeChild.data.crossLinks) {
              if (dRoot.data.key.startsWith(link[0])) {
                label = link[1]
                break
              }
            }
            if (label) break
          }
          crossLinks.push({source: child, target: dRoot, id: `${child.id}-${dRoot.id}`, label: label})
        }
      })
      dRoot._children = dRoot.children
      node.subGraph = {root: dRoot, crossLinks: crossLinks}
      subGraphNodes.push(...dRoot.descendants())
    }

    // TODO using hierarchy is not required
    const root = d3.hierarchy(data)
    nodes = root.descendants()

    const subGraphNodes = []

    const crossLinkMap = {}

    nodes.forEach((d, i) => {
      // initialize all nodes at the center
      d.x = width / 2
      d.y = width / 2
      d.id = i
      d._children = d.children
      if (d.depth) d.children = null
      if (d.id) if (!crossLinkMap[d.data.key]) crossLinkMap[d.data.key] = d
      if (d.data.type && d.data.type.startsWith('task')) addSubGraph(d)
    })
    addMarkers(svg, nodes, circleRadius)

    addMarkers(svg, subGraphNodes)

    // add zoom
    const zoomBehaviors = d3.zoom()
    .on('zoom', () => {
      svgGroup.attr('transform', d3.event.transform)
    })

    svg.call(zoomBehaviors)

    const simulation = d3.forceSimulation()
      .force('crosslink', d3.forceLink().id(d => d.id).distance(linkDistance).strength(1.0))
      .force('link', d3.forceLink().id(d => d.id).distance(linkDistance).strength(0.5))
      .force('charge', d3.forceManyBody().strength(-50.0))
      .velocityDecay(0.2)
      .alphaDecay(0.2)
      .force('x', d3.forceX(d => d.data.fixX || 0).strength(d => d.data.fixX === undefined ? 0 : 3))
      .force('y', d3.forceY(d => d.data.fixY || 0).strength(d => d.data.fixY === undefined ? 0 : 3))
      .force('collision', d3.forceCollide().radius(circleRadius * 2))

    function drag(simulation) {
      function dragStart(d) {
        if (!d3.event.active) simulation.alphaTarget(0.7).restart()
        d.fx = d.x
        d.fy = d.y
      }

      function dragged(d) {
        d.fx = d3.event.x
        d.fy = d3.event.y
      }

      function dragEnd(d) {
        if (!d3.event.active) simulation.alphaTarget(0)
        d.fx = null
        d.fy = null
      }

      return d3.drag()
        .on('start', dragStart)
        .on('drag', dragged)
        .on('end', dragEnd)
    }

    let node, link, crossLink

    function ticked() {
      link
        .attr('x1', d => d.source.x)
        .attr('y1', d => d.source.y)
        .attr('x2', d => d.target.x)
        .attr('y2', d => d.target.y)

      crossLink
        .attr('x1', d => d.source.x)
        .attr('y1', d => d.source.y)
        .attr('x2', d => d.target.x)
        .attr('y2', d => d.target.y)

      node
        .attr('transform', d => `translate(${d.x},${d.y})`)
    }

    let nodesHistory = []
    let currentNodesIndex = 0
    const nodesMap = {}
    nodes.forEach(node => {
      nodesMap[node.id] = node
    })

    const pushNodes = (nodes) => {
      nodesHistory.push(nodes.map(node => node.id))
      currentNodesIndex = currentNodesIndex + 1
    }

    const rootNodes = root.leaves()
    pushNodes(rootNodes)

    d3.select('#backbutton')
      .on('click', () => {
        currentNodesIndex = currentNodesIndex - 1
        currentNodesIndex = Math.max(currentNodesIndex, 0)
        update(nodesHistory[currentNodesIndex].map(id => nodesMap[id]))
        simulation.restart()
      })

    d3.select('#forwardbutton')
      .on('click', () => {
        currentNodesIndex = currentNodesIndex + 1
        currentNodesIndex = Math.min(currentNodesIndex, nodesHistory.length - 1)
        update(nodesHistory[currentNodesIndex].map(id => nodesMap[id]))
        simulation.restart()
      })

    d3.select('#resetbutton')
      .on('click', () => {
      currentNodesIndex = 0
      nodesHistory = []
      pushNodes(rootNodes)
      update(rootNodes)
      simulation.restart()
    })

    function update(nodes) {
      const links = []
      const crossLinks = []
      nodes.forEach(node => {
        const nodeLinks = node.data.crossLinks || []
        nodeLinks.forEach(link => {
          const nodeTarget = crossLinkMap[link[0]]
          if (node && nodes.includes(nodeTarget)) {
            crossLinks.push({
              source: node.id,
              target: nodeTarget.id,
              label: link[1],
              id: `${node.id}-${nodeTarget.id}`,
              dash: link[2],
              marker: link[3] === undefined
            })
          }
        })
      })

      //  update nodes
      node = gNode.selectAll('g')
        .data(nodes, d => d.id)

      node.exit().remove()

      const handleMouseOverText = d => {
        d3.select(`#text-${d.id}`).style('fill', red[800])
        const text = d.data.type.includes('workflow') ? 'overview page' : 'archive section'
        if (d.data.entryId) {
          setTooltipContent(`Go to ${text} for entry ${d.data.entryId}`)
        }
      }

      const handleMouseOutText = d => {
        setTooltipContent('')
        d3.select(`#text-${d.id}`).style('fill', grey[800])
      }

      const handleClickText = d => {
        if (d.id !== 0) {
          let path = `entry/id/${d.data.entryId}`
          const sectionPath = d.data.path ? d.data.path.replace(/\/(?=\d)/g, ':') : null
          path = d.data.type.startsWith('workflow') ? path : sectionPath ? `${path}/data${sectionPath}` : path
          const url = getUrl(path)
          history.push(url)
        }
      }

      const handleMouseOverCircle = d => {
        d3.select(`#circle-${d.id}`).style('stroke-opacity', 1)
        const children = d.children || []
        children.forEach(child => {
          d3.select(`#link-${d.id}-${child.id}`).style('stroke-opacity', 1)
          d3.select(`#circle-${child.id}`).style('stroke-opacity', 1)
        })
      }

      const handleMouseOutCircle = d => {
        d3.select(`#circle-${d.id}`).style('stroke-opacity', 0.5)
        const children = d.children || []
        children.forEach(child => {
          d3.select(`#link-${d.id}-${child.id}`).style('stroke-opacity', 0.5)
          d3.select(`#circle-${child.id}`).style('stroke-opacity', 0.5)
        })
      }

      const handleMouseOverCrossLink = d => {
        d3.select(`#crosslink-${d.id}`).style('stroke-opacity', 1.0)
        d3.select(`.marker-${d.id}`).style('fill-opacity', 1.0)
        d3.select(`#circle-${d.source.id}`).style('stroke', red[500]).style('stroke-opacity', 1.0)
        d3.select(`#circle-${d.target.id}`).style('stroke', red[500]).style('stroke-opacity', 1.0)
        setTooltipContent(d.label)
      }
      const handleMouseOutCrossLink = d => {
        d3.select(`#crosslink-${d.id}`).style('stroke-opacity', 0.5)
        d3.select(`.marker-${d.id}`).style('fill-opacity', 0.5)
        d3.select(`#circle-${d.source.id}`).style('stroke', grey[500]).style('stroke-opacity', 0.5)
        d3.select(`#circle-${d.target.id}`).style('stroke', grey[500]).style('stroke-opacity', 0.5)
        setTooltipContent('')
      }

      const nodeEnter = node.enter().append('g')
        .attr('id', d => `node-${d.id}`)
        .on('click', d => {
          if (d3.event.defaultPrevented) return

          if (d.data.type.startsWith('task') && d.subGraph) {
            const k = width / (circleRadius * 4)
            const dNodes = d.subGraph.root.descendants()

            gNode.selectAll('g').selectAll('.subgraph')
              .attr('visibility', g => g.id === d.id ? 'visible' : 'hidden')

            const node = d3.select(`#node-${d.id}`)
              .raise()

            node.select('.subgraph')
              .attr('visibility', 'visible')

            const dNode = node.select('.subgraph')
              .raise()

            const subNode = dNode.selectAll('circle')
              .data(d.subGraph.root.children ? dNodes : [], d => d.id)

            const subNodeEnter = subNode.enter()
              .append('circle')
              .attr('class', 'circle')
              .attr('id', d => `circle-${d.id}`)
              .attr('r', di => di.r * k)
              .attr('fill', d => d.depth ? nomadPrimaryColor.light : nomadSecondaryColor.dark)
              .attr('transform', di => `translate(${(di.x - dNodes[0].x) * k},${(di.y - dNodes[0].y) * k})`)
              .on('mouseover', handleMouseOverCircle)
              .on('mouseout', handleMouseOutCircle)

            subNode.exit().remove()

            subNodeEnter.merge(subNode)

            const dCrossLinks = d.subGraph.root.children ? d.subGraph.crossLinks : []
            const subCrossLink = dNode.select('.crosslink').raise().selectAll('path')
              .data(dCrossLinks, d => d.id)

            const subCrossLinkEnter = subCrossLink.enter()
              .append('path')
              .attr('class', 'crosslink')
              .attr('id', d => `crosslink-${d.id}`)
              .attr('marker-end', d => `url(#${d.id})`)
              .attr('d', d => {
                  const fcrossLink = (source, target) => {
                    const path = d3.path()
                    path.moveTo(source.x, source.y)
                    path.lineTo(target.x, target.y)
                    return path
                  }
                  const translate = (node) => {
                    return {x: (node.x - dNodes[0].x) * k, y: (node.y - dNodes[0].y) * k}
                  }
                  return fcrossLink(translate(d.source), translate(d.target))
              })
              .on('mouseover', handleMouseOverCrossLink)
              .on('mouseout', handleMouseOutCrossLink)

            subCrossLink.exit().remove()

            subCrossLinkEnter.merge(subCrossLink)

            const subLabel = dNode.selectAll('text')
              .data(d.subGraph.root.children ? dNodes : [], d => d.id)

            const subLabelEnter = subLabel.enter()
              .append('text')
              .attr('class', 'text')
              .attr('id', d => `text-${d.id}`)
              .style('font-size', d => d.parent ? 12 : 18)
              .attr('y', d => d.parent ? -10 : -d.r * k * 0.8)
              .attr('text-anchor', 'middle')
              .text(d => d.data.name)
              .attr('transform', di => `translate(${(di.x - dNodes[0].x) * k},${(di.y - dNodes[0].y) * k})`)
              .on('mouseover', handleMouseOverText)
              .on('mouseout', handleMouseOutText)
              .on('click', handleClickText)

            subLabel.exit().remove()

            subLabelEnter.merge(subLabel)

            d.subGraph.root.children = d.subGraph.root.children ? null : d.subGraph.root._children
          } else if (d._children) {
            if (d.children) {
              d.children = null
            } else {
              const tasks = d._children.filter(child => {
                return child.data.type && (child.data.type.startsWith('workflow') || child.data.type.startsWith('task'))
              })
              if (tasks.length > 0) nodes = nodes.filter(node => node.id !== d.id)
              const keys = nodes.map(node => node.data.key)
              d._children.forEach(child => {
                if (!keys.includes(child.data.key)) nodes.push(child)
              })
              d.children = d._children
            }
            // save snapshots of nodes
            pushNodes(nodes)
            update(nodes)
          }
        })
        .call(drag(simulation))

      const subGraph = nodeEnter.append('g')
        .attr('class', 'subgraph')

      subGraph.append('g')
        .attr('class', 'crosslink')

      nodeEnter.append('circle')
        .attr('class', 'circle')
        .attr('id', d => `circle-${d.id}`)
        .attr('r', circleRadius)
        .attr('fill', d => {
          if (d.data.type.startsWith('workflow')) return nomadPrimaryColor.dark
          if (d.data.type.startsWith('task')) return nomadSecondaryColor.dark
          return nomadPrimaryColor.light
        })
        .on('mouseover', handleMouseOverCircle)
        .on('mouseout', handleMouseOutCircle)

      nodeEnter.append('text')
        .attr('class', 'text')
        .attr('id', d => `text-${d.id}`)
        .text(d => d.data.name)
        .style('font-size', d => d.depth ? 12 : 18)
        .attr('text-anchor', 'middle')
        .attr('y', -circleRadius - 5)
        .on('click', handleClickText)
        .on('mouseover', handleMouseOverText)
        .on('mouseout', handleMouseOutText)

      // update the text
      node.selectAll('text')
        .text(d => d.data.name)

      node = nodeEnter.merge(node)

      link = gLink.selectAll('line')
        .data(links, d => d.target.id)

      link.exit().remove()

      const linkEnter = link.enter().append('line')
        .attr('id', d => `link-${d.source.id}-${d.target.id}`)

      link = linkEnter.merge(link)

      crossLink = gCrossLink.selectAll('line')
        .data(crossLinks, d => d.id)

      crossLink.exit().remove()

      const crossLinkEnter = crossLink.enter().append('line')
        .attr('class', 'crosslink')
        .style('stroke-dasharray', d => d.dash)
        .attr('marker-end', d => d.marker ? `url(#${d.id})` : null)
        .attr('id', d => `crosslink-${d.id}`)
        .on('mouseover', handleMouseOverCrossLink)
        .on('mouseout', handleMouseOutCrossLink)

      crossLink = crossLinkEnter.merge(crossLink)

      simulation
        .nodes(nodes)
        .on('tick', ticked)

      // disabled links links
      // simulation
      //   .force('link')
      //   .links(links)

      simulation
        .force('crosslink')
        .links(crossLinks)
    }

    update(rootNodes)
  }, [data, svgRef, history, finalLayout, setTooltipContent])
  return <svg className={classes.root} ref={svgRef}></svg>
})

ForceDirected.propTypes = {
  data: PropTypes.object.isRequired,
  layout: PropTypes.object,
  setTooltipContent: PropTypes.any
}

const WorkflowCard = React.memo(({archive}) => {
  const [data, setData] = useState()
  const {api} = useApi()
  const [tooltipContent, setTooltipContent] = useState('')
  const [tooltipPosition, setTooltipPosition] = useState({x: undefined, y: undefined})

  useEffect(() => {
    const crossLinks = {}

    const addCrossLink = (key, link) => {
      if (key in crossLinks) {
        crossLinks[key].push(link)
      } else {
        crossLinks[key] = [link]
      }
    }

    const createHierarchy = async function(section, archive, name, type, index, reference) {
      const baseUrl = createEntryUrl(apiBase, archive?.metadata?.upload_id, archive?.metadata?.entry_id)

      if (typeof section === 'string') {
        const match = section.match('.*/(\\w+)/(\\d+)$')
        if (match) {
          type = type || match[1]
          index = match.length > 2 ? parseInt(match[2]) : 0
        } else {
          type = type || section.split('/').pop()
        }
        reference = section
      }

      // resolve section from reference path
      const resolved = await (async (path) => {
        if (!(typeof path === 'string')) {
          return [section, archive, path, baseUrl]
        }
        if (path.startsWith('#')) path = path.slice(1)
        if (!path.includes('#/')) {
          try {
            return [resolveInternalRef(path, archive), archive, path, baseUrl]
          } catch (error) {
            return [null, archive, path, baseUrl]
          }
        }
        const url = resolveNomadUrl(path, baseUrl)
        const query = {'workflow2': '*', 'metadata': '*'}
        try {
          const response = await api.post(`/entries/${url.entryId}/archive/query`, {required: query})
          let sectionPath = path.split('#').pop()
          if (!sectionPath.startsWith('/')) sectionPath = `/${sectionPath}`
          const archive = response.data.archive
          const section = type !== 'section' ? resolveInternalRef(sectionPath, archive) : sectionPath
          const baseUrl = createEntryUrl(apiBase, archive?.metadata?.upload_id, archive?.metadata?.entry_id)
          return [section, archive, sectionPath, baseUrl]
        } catch (error) {
          console.error(`Cannot resolve entry ${url.entryId}: ${error}`)
          return [null, archive, path, baseUrl]
        }
      })(reference)

      // get data for the current section
      const sectionData = await (async (section, archive, path, baseUrl) => {
        const sectionKey = [baseUrl, path].join('#')
        const children = []
        const maxNodes = 6
        const start = Math.floor(maxNodes / 2)
        if (typeof section === 'string' || !section) {
          return {
            name: name,
            key: sectionKey,
            type: type,
            size: 20,
            children: null,
            path: path.split('#').pop(),
            entryId: archive?.metadata?.entry_id,
            index: index
          }
        }
        if (Array.isArray(section)) {
          for (const [index, ref] of section.entries()) {
            const child = await createHierarchy(ref, archive, ref.name, type, index, `${path}/${index}`)
            children.push(child)
          }
        } else {
          const taskInputs = []
          const taskOutputs = []
          if (section.section) {
            const child = await createHierarchy(section.section, archive, section.name, 'section')
            children.push(child)
          }
          if (section.inputs) {
            let parent
            const end = section.inputs.length - start
            for (const [index, ref] of section.inputs.entries()) {
              if (index < start || index >= end || start + 1 === end) {
                parent = await createHierarchy(ref, archive, ref.name || 'input', 'inputs', index, `${path}/inputs/${index}`)
                if (!parent) continue
                const child = parent.children ? parent.children[0] : null
                if (child) {
                  addCrossLink(child.key, [sectionKey, ref.name || ''])
                  child.intent = 'input'
                  taskInputs.push(child)
                }
              }
            }
            if (section.inputs.length > maxNodes + 1) {
              const otherInputs = {
                name: `Inputs ${start + 1} - ${end} not shown`,
                type: taskInputs[start - 1].type,
                path: taskInputs[start - 1].path,
                entryId: taskInputs[start - 1].entryId,
                key: `${taskInputs[start - 1].key}.invisible.input`,
                intent: 'input',
                size: 20
              }
              taskInputs.splice(start, 0, otherInputs)
              addCrossLink(taskInputs[start - 1].key, [otherInputs.key, '', null, false])
              addCrossLink(taskInputs[start + 1].key, [otherInputs.key, '', null, false])
            }
            children.push(...taskInputs)
          }
          if (section.outputs) {
            const end = section.outputs.length - start
            for (const [index, ref] of section.outputs.entries()) {
              if (index < start || index >= end || start + 1 === end) {
                const parent = await createHierarchy(ref, archive, ref.name || 'output', 'outputs', index, `${path}/outputs/${index}`)
                if (!parent) continue
                const child = parent.children ? parent.children[0] : null
                if (child) {
                  addCrossLink(sectionKey, [child.key, ref.name || ''])
                  child.intent = 'output'
                  taskOutputs.push(child)
                }
              }
            }
            if (section.outputs.length > maxNodes + 1) {
              const otherOutputs = {
                name: `Outputs ${start + 1} - ${end} not shown`,
                type: taskOutputs[start - 1].type,
                path: taskOutputs[start - 1].path,
                entryId: taskOutputs[start - 1].entryId,
                key: `${taskOutputs[start - 1].key}.invisible.output`,
                intent: 'output',
                size: 20
              }
              taskOutputs.splice(start, 0, otherOutputs)
              addCrossLink(taskOutputs[start - 1].key, [otherOutputs.key, '', null, false])
              addCrossLink(taskOutputs[start + 1].key, [otherOutputs.key, '', null, false])
            }
            children.push(...taskOutputs)
          }
          if (section.task) {
            const task = await createHierarchy(section.task, archive)
            task.name = section.name || task.name
            task.children = task.children || []
            const taskKeys = task.children.map(child => child.key)
            taskInputs.forEach(child => {
              if (!taskKeys.includes(child.key)) task.children.push(child)
            })
            taskOutputs.forEach(child => {
              if (!taskKeys.includes(child.key)) task.children.push(child)
            })
            children.push(task)
            return task
          }
          if (section.tasks) {
            type = 'workflow'
            // resolve only several first and last tasks
            const end = section.tasks.length - start
            const sectionChildren = []
            for (const [index, ref] of section.tasks.entries()) {
              if (index < start || index >= end || start + 1 === end) {
                const task = await createHierarchy(ref, archive, ref.name || 'task', 'tasks', index, `${path}/tasks/${index}`)
                if (task) sectionChildren.push(task)
              }
            }
            if (section.tasks.length > maxNodes + 1) {
              const otherTasks = {
                name: `Tasks ${start + 1} - ${end} not shown`,
                type: sectionChildren[start - 1].type,
                path: sectionChildren[start - 1].path,
                entryId: sectionChildren[start - 1].entryId,
                key: `${sectionChildren[start - 1].key}.invisible.task`,
                size: 20
              }
              sectionChildren.splice(start, 0, otherTasks)
              addCrossLink(sectionChildren[start - 1].key, [otherTasks.key, '', null, false])
              addCrossLink(sectionChildren[start + 1].key, [otherTasks.key, '', null, false])
            }
            children.push(...sectionChildren)
          }
        }

        return {
          name: section.name || name || type,
          key: sectionKey,
          type: type,
          size: 20,
          children: children.length === 0 ? null : children,
          path: path,
          entryId: archive?.metadata?.entry_id,
          index: index
        }
      })(resolved[0], resolved[1], resolved[2], resolved[3])

      return sectionData
    }

    const workflow = archive?.workflow2
    if (workflow) {
      createHierarchy('/workflow2', archive, 'Workflow').then(result => {
        // save the cross links to each section
        function getLinks(section) {
          section.crossLinks = crossLinks[section.key]
          if (section.children) {
            section.children.forEach(child => {
            getLinks(child)
            })
          }
        }
        getLinks(result)
        setData(result)
      })
    }
  }, [archive, api])

  const graph = useMemo(() => {
    if (!data) {
      return ''
    }

    return <ForceDirected data={data} setTooltipContent={setTooltipContent}></ForceDirected>
  }, [data])

  const actions = <div>
    <IconButton id='backbutton'>
      <Tooltip title="Back">
        <Undo />
      </Tooltip>
    </IconButton>
    <IconButton id='forwardbutton'>
      <Tooltip title="Forward">
        <Redo />
      </Tooltip>
    </IconButton>
    <IconButton id='resetbutton'>
      <Tooltip title="Reset">
        <Replay />
      </Tooltip>
    </IconButton>
  </div>

  return graph && <PropertyCard title='Workflow Graph' action={actions}>
    <PropertyGrid>
      <Tooltip title={tooltipContent}
        onMouseMove={event => setTooltipPosition({x: event.pageX, y: event.pageY})}
        enterNextDelay={700}
        leaveDelay={700}
        arrow
        PopperProps={
          {anchorEl: {
            clientHeight: 0,
            clientWidth: 0,
            getBoundingClientRect: () => ({
              top: tooltipPosition.y,
              left: tooltipPosition.x,
              right: tooltipPosition.x,
              bottom: tooltipPosition.y,
              width: 0,
              height: 0
            })
          }}
        }
        >
          <div>
            {graph}
          </div>
        </Tooltip>
    </PropertyGrid>
  </PropertyCard>
})

WorkflowCard.propTypes = {
  index: PropTypes.object,
  properties: PropTypes.object,
  archive: PropTypes.object
}

export default WorkflowCard
