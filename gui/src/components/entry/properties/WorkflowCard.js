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
import * as d3 from 'd3'
import { makeStyles, Tooltip, IconButton, TextField, FormControl } from '@material-ui/core'
import grey from '@material-ui/core/colors/grey'
import blueGrey from '@material-ui/core/colors/blueGrey'
import red from '@material-ui/core/colors/red'
import { Replay, Undo, Label, LabelOff, PlayArrowSharp, StopSharp, Clear } from '@material-ui/icons'
import { useHistory } from 'react-router-dom'
import { isPlainObject } from 'lodash'
import { PropertyCard, PropertyGrid } from './PropertyCard'
import { resolveNomadUrl, resolveInternalRef, createEntryUrl } from '../../../utils'
import { useApi } from '../../api'
import { getUrl } from '../../nav/Routes'
import { nomadFontFamily, apiBase } from '../../../config'
import { useAsyncError } from '../../../hooks'

const useWorkflowGraphStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    minWidth: 1000,
    '& .link': {
      fill: 'none',
      strokeOpacity: 0.5,
      fillOpacity: 0.5,
      strokeWidth: 2.5
    },
    '& .text': {
      fontFamily: nomadFontFamily,
      fontSize: 12,
      dy: '0.35em',
      textAnchor: 'middle'
    },
    '& .icon': {
      strokeWidth: 2,
      strokeOpacity: 0.5
    }
  }
}))

let archives = {}

const resolveSection = async (source, query) => {
  if (isPlainObject(source.section)) {
    // already resolved
    return source
  }

  let path = source.path
  if (typeof path !== 'string') path = source.section || '/'
  const pathSegments = path.split('#').filter(p => p)

  let archive = source.archive
  let baseUrl = createEntryUrl(apiBase, archive?.metadata?.upload_id, archive?.metadata?.entry_id)
  let section
  if (pathSegments.length === 1) {
    // internal reference
    path = pathSegments[0]
    try {
      section = source.section ? source.path : resolveInternalRef(path, archive)
    } catch (error) {
      console.error(`Cannot resolve section ${path}`)
      return
    }
  } else {
    // external reference
    // temporary fix for path following the format /entries/entry_id/archive
    // TODO extend resolveNomadUrl for other reference formats
    let entryId
    const match = path.match('/entries/(.+?)/archive')
    if (match) {
      entryId = match[1]
    } else {
      const url = resolveNomadUrl(path, baseUrl)
      entryId = url.entryId
    }
    const {api, required} = query
    try {
      archive = archives[entryId]
      if (!archive) {
        const response = await api.post(`/entries/${entryId}/archive/query`, {required: required})
        archive = response.data.archive
        archives[entryId] = archive
      }
      path = pathSegments[1]
      if (!path.startsWith('/')) path = `/${path}`
      section = source.section ? source.path : resolveInternalRef(path, archive)
    } catch (error) {
      console.error(`Cannot resolve entry ${entryId}: ${error}`)
      return
    }
  }

  if (!section) return

  baseUrl = createEntryUrl(apiBase, archive?.metadata?.upload_id, archive?.metadata?.entry_id)
  const match = path.match('.*/([^/]+)/(\\d+)$')

  const sectionKeys = query.sectionKeys || []
  const nChildren = (section) => (sectionKeys.map(key => (section[key] || []).length)).reduce((a, b) => a + b)
  if (section.task) {
    const task = await resolveSection({path: section.task, archive: archive}, query)
    if (!task) return
    task.section.inputs = section.inputs && section.inputs.length ? section.inputs : task.section.inputs
    task.section.outputs = section.outputs && section.outputs.length ? section.outputs : task.section.outputs
    task.nChildren = nChildren(task.section)
    return task
  }

  return {
    name: section.name,
    section: section,
    sectionType: match ? match[1] : path.split('/').pop(),
    path: path,
    archive: source.section ? null : archive,
    url: [baseUrl, path].join('#'),
    entryId: archive?.metadata?.entry_id,
    color: section.color,
    nChildren: nChildren(section)
  }
}

const range = (start, stop) => Array.from(Array(stop - start).fill(start).map((n, i) => n + i))

const getNodes = async (source, query) => {
  let {resolveIndices, maxNodes} = query
  const {sectionKeys} = query

  let nodes = source.nodes || []
  if (nodes.length && !resolveIndices) return nodes

  if (!maxNodes) maxNodes = 7

  const resolved = await resolveSection(source, query)
  if (!resolved) return nodes
  const parent = resolved.section
  if (!parent) return nodes

  nodes = []
  const archive = resolved.archive
  for (const key of sectionKeys) {
    let children = parent[key] || []
    const nChildren = children.length
    if (!Array.isArray(children)) children = [children]

    let sectionIndices = range(0, nChildren)
    if (resolveIndices && resolveIndices[key]) {
      sectionIndices = resolveIndices[key] || sectionIndices
    } else if (maxNodes < nChildren) {
      const mid = Math.floor(maxNodes / 2)
      // show first and last few nodes
      sectionIndices = [...range(0, mid), ...range(nChildren - mid, nChildren)]
    }

    for (const index of sectionIndices) {
      const child = children[index]
      if (!child) continue

      let path = child
      if (!child.section) path = `${source.path}/${key}/${index}`
      const section = await resolveSection({section: child.section, path: path, archive: archive}, query)
      if (!section) continue

      section.name = child.name
      section.type = key
      section.index = index
      section.parent = source
      section.total = nChildren
      nodes.push(section)
    }
  }

  return nodes
}

const getLinks = async (source, query) => {
  if (source.links) return source.links

  const nodes = source.nodes
  if (!source.nodes) return

  const links = []
  for (const node of nodes) {
    if (node.type === 'tasks') {
      node.nodes = await getNodes(node, query)
    }
  }

  const isLinked = (source, target) => {
    if (source.url === target.url) return false

    const outputs = []
    if (source.type === 'tasks' && source.nodes) {
      outputs.push(...source.nodes.filter(node => node.type === 'outputs').map(node => node.url))
    } else {
      outputs.push(source.url)
    }

    const inputs = []
    if (target.type === 'tasks' && target.nodes) {
      inputs.push(...target.nodes.filter(node => node.type && node.type.startsWith('inputs')).map(node => node.url))
    } else {
      inputs.push(target.url)
    }

    let linked = false
    for (const output of outputs) {
      if (!output) continue
      if (inputs.includes(output)) {
        linked = true
        break
      }
    }
    return linked
  }

  // links from inputs to source
  const inputs = nodes.filter(node => node.type && node.type.startsWith('inputs') && node.url)
  inputs.forEach(input => {
    links.push({source: input, target: source, label: 'Input'})
  })
  // links from source to outputs
  const outputs = nodes.filter(node => node.type === 'outputs' && node.url)
  outputs.forEach(output => {
    links.push({source: source, target: output, label: 'Output'})
  })

  // source and target are linked if any of the outputs of source or source itself is any of the
  // inputs of target or target itself
  for (const source of ['inputs', 'tasks']) {
    for (const target of ['outputs', 'tasks']) {
      const sourceNodes = nodes.filter(node => node.type && node.type.startsWith(source) && node.url)
      const targetNodes = nodes.filter(node => node.type === target && node.url)
      sourceNodes.forEach(sourceNode => {
        targetNodes.forEach(targetNode => {
          if (isLinked(sourceNode, targetNode)) {
            let label = ''
            if (source === 'inputs') {
              label = 'Input'
            } else if (target === 'outputs') {
              label = 'Output'
            } else {
              label = 'Click to see how tasks are linked'
            }
            links.push({source: sourceNode, target: targetNode, label: label})
          }
        })
      })
    }
  }

  return links
}

const Graph = React.memo(({
  source,
  query,
  layout,
  setTooltipContent,
  setCurrentNode,
  setShowLegend,
  setEnableForce
  }) => {
  const classes = useWorkflowGraphStyles()
  const svgRef = useRef()
  const history = useHistory()
  const asyncError = useAsyncError()
  const finalLayout = useMemo(() => {
    const defaultLayout = {
      width: 700,
      margin: {top: 60, bottom: 60, left: 40, right: 40},
      circleRadius: 20,
      markerWidth: 4,
      scaling: 0.35,
      color: {
        text: grey[800],
        link: red[800],
        outline: blueGrey[800],
        workflow: '#192E86',
        task: '#00AC7C',
        input: '#A59FFF',
        output: '#005A35'
      },
      shape: {
        input: 'circle',
        workflow: 'rect',
        task: 'rect',
        output: 'circle'
      }
    }
    return {...defaultLayout, ...layout}
  }, [layout])
  archives = {}

  useEffect(() => {
    const {width, markerWidth, scaling, color} = finalLayout
    const nodeShape = finalLayout.shape
    const whRatio = 1.6
    const legendSize = 12
    const height = width / whRatio

    const svg = d3.select(svgRef.current)

    const inOutColor = d3.interpolateRgb(color.input, color.output)(0.5)

    let nodes, node, line
    let id = 0
    let view
    let focus
    let previousNode = 'root'
    let root
    const dquery = {...query, resolveIndices: {}}
    let currentNode = root
    let hasError = false
    let enableForce = false
    let dragged = false
    let zoomTransform = d3.zoomIdentity

    const isWorkflow = (d) => {
      if (d.sectionType && d.sectionType.startsWith('workflow')) return true
      const tasks = (d.nodes || []).filter(n => n.type === 'tasks')
      return tasks.length > 0
    }

    const nodeColor = (d) => {
      if (d.color) return d.color
      if (d.type === 'link') return '#ffffff'
      if (isWorkflow(d)) return color.workflow
      if (d.type === 'tasks') return color.task
      if (d.type === 'inputs-outputs') return inOutColor
      return d.type === 'inputs' ? color.input : color.output
    }

    const trimName = (name) => name && name.length > 25 ? `${name.substring(0, 22)}...` : name

    svg.attr('width', width)
      .attr('height', height)

    svg.append('g')
      .attr('class', 'svgGroup')

    const svgGroup = svg.select('.svgGroup')

    svg.append('defs')
      .attr('class', 'defs')

    const addLinkMarkers = (links) => {
      svg.select('.defs')
        .selectAll('marker')
        .exit().remove()
        .data(links.map(link => link.id))
        .enter()
        .append('marker')
        .attr('class', d => `marker-${d}`)
        .attr('id', String)
        .attr('viewBox', '0 -5 10 10')
        .attr('refX', 0)
        .attr('refY', 0)
        .attr('markerWidth', markerWidth)
        .attr('markerHeight', markerWidth)
        .attr('xoverflow', 'visible')
        .attr('orient', 'auto')
        .style('fill', color.link)
        .attr('fill-opacity', 0.5)
        .append('path')
        .attr('d', 'M0,-5L10,0L0,5')
    }

    svg.append('g')
      .attr('class', 'legend')
      .attr('visibility', 'visible')

    const legend = svg.select('.legend')

    const addLegend = (label, index) => {
      const gLegend = legend.append('g')
        .attr('cursor', 'pointer')
        .attr('pointer-events', 'all')
      const shape = nodeShape[label]
      const icon = gLegend.append(shape)
      const dx = width / 8
      const x = width / 2 + dx * index - dx * 1.5
      const y = height - legendSize * 2

      if (shape === 'rect') {
        icon.attr('width', legendSize * whRatio)
          .attr('height', legendSize)
          .attr('rx', legendSize * 0.1)
          .attr('x', x - legendSize * whRatio / 2)
          .attr('y', y - legendSize / 2)
      } else {
        icon.attr('r', legendSize / 2)
          .attr('cx', x)
          .attr('cy', y)
        }
      icon
        .attr('fill', nodeColor({type: label + 's', sectionType: label + 's'}))

      gLegend.append('text')
        .attr('class', 'text')
        .text(label)
        .style('text-anchor', 'middle')
        .attr('x', x)
        .attr('y', y + legendSize * 1.2)
        .style('alignment-baseline', 'middle')

      gLegend
        .on('mouseover', () => {
          let tooltip = ''
          if (label === 'input') {
            tooltip = <p>
              Input to a task or workflow.
            </p>
          } else if (label === 'output') {
            tooltip = <p>
              Output from a task or workflow.
            </p>
          } else if (label === 'workflow') {
            tooltip = <p>
              Task containing further sub-tasks.
            </p>
          } else if (label === 'task') {
            tooltip = <p>
              Elementary task with inputs and outputs.
            </p>
          }
          setTooltipContent(tooltip)
        })
        .on('mouseout', () => {
          setTooltipContent('')
        })
    }

    // add legend
    legend.append('path')
      .attr('d', () => {
        const line = d3.line().x(d => d[0]).y(d => d[1])
        return line([
          [width / 4 + 70, height - legendSize * 3],
          [3 * width / 4, height - legendSize * 3],
          [3 * width / 4, height],
          [width / 4, height],
          [width / 4, height - legendSize * 3],
          [width / 4 + 10, height - legendSize * 3]
        ])
      })
      .attr('stroke', grey[900])
      .attr('stroke-width', 0.5)
      .attr('fill', 'none')
    legend.append('text')
      .text('Legend')
      .attr('class', 'text')
      .attr('font-weight', 'bold')
      .attr('x', width / 4 + 40)
      .attr('y', height - legendSize * 2.8)

    const legendLabels = Object.keys(nodeShape)
    legendLabels.forEach((label, index) => addLegend(label, index))

    // add zoom
    const zoomBehaviors = d3.zoom()
      .on('zoom', () => {
        zoomTransform = d3.event.transform
        svgGroup.attr('transform', zoomTransform)
      })

    svg.call(zoomBehaviors)

    // set zoom factor to <1 inorder to see the inputs and outputs
    const zoomF = (r) => width * scaling / r

    // add force directed behavior for tasks inside box
    const simulation = d3.forceSimulation()
      .force('link', d3.forceLink().id(d => d.id))
      .force('charge', d3.forceManyBody())
      .velocityDecay(0.2)
      .alphaDecay(0.2)

    // add drag
    const dragBehaviors = d3.drag()
      .on('drag', d => {
        if (!d.parent || d.id === nodes[0].id) return
        const event = d3.event.sourceEvent
        const k = zoomF(view[2])
        const x = ((event.offsetX - zoomTransform.x) / zoomTransform.k - width / 2) / k + view[0]
        const y = ((event.offsetY - zoomTransform.y) / zoomTransform.k - height / 2) / k + view[1]
        setTooltipContent('')
        if (enableForce) {
          d.fx = x
          d.fy = y
        } else {
          d.x = x
          d.y = y
          zoomTo(view)
        }
      })
      .on('start', d => {
        if (!enableForce) return
        simulation.alphaTarget(0.7).restart()
      })
      .on('end', d => {
        dragged = true
        if (!enableForce) return
        simulation.alphaTarget(0)
        d.fx = null
        d.fy = null
      })

    const gNode = svgGroup.append('g')
      .attr('class', 'node')
      .attr('cursor', 'pointer')
      .attr('pointer-events', 'all')

    const gLink = svgGroup.append('g')
      .attr('class', 'link')
      .attr('cursor', 'default')
      .attr('pointer-events', 'all')
      .attr('stroke', color.link)

    const setNodesPosition = (root) => {
      // layout root and its nodes horizontally from inputs (left) root and tasks (middle)
      // and outputs (right)
      if (!root.nodes || root.children) return

      const circleRadius = root.r / 12
      const dx = circleRadius * 4
      const fx = 1.4
      const tasks = root.nodes.filter(node => node.type === 'tasks')
      root.children = tasks
      // tasks are arranged in a circle inside the parent workflow
      tasks.forEach((node, index) => {
        const theta = (2 * Math.PI * index / tasks.length) - Math.PI / 4
        node.r = Math.sqrt(node.size || 1) * circleRadius * 1.5
        const r = root.r - node.r
        node.x = -r * 0.95 * Math.cos(theta) * whRatio + root.x
        node.y = -r * (theta < 0 || theta > Math.PI ? 0.95 : 0.85) * Math.sin(theta) + root.y
      })

      // inputs left
      const inputs = root.nodes.filter(node => node.type === 'inputs')
      let offsetX = (inputs.length - 1) * dx / 2
      inputs.forEach((node, index) => {
        node.y = index * dx - offsetX + root.y
        node.x = -root.r * fx * whRatio + root.x
        node.r = circleRadius
      })

      // inputs-outputs top
      const inouts = root.nodes.filter(node => node.type === 'inputs-outputs')
      offsetX = (inouts.length - 1) * dx
      inouts.forEach((node, index) => {
        node.x = index * dx * 2 - offsetX + root.x
        node.y = -root.r * fx + root.y
        node.r = circleRadius
      })

      // outputs right
      const outputs = root.nodes.filter(node => node.type === 'outputs')
      offsetX = (outputs.length - 1) * dx / 2
      outputs.forEach((node, index) => {
        node.y = index * dx - offsetX + root.y
        node.x = root.r * whRatio * fx + root.x
        node.r = circleRadius
      })
      root.nodes = [...inputs, ...inouts, ...root.children, ...outputs]
    }

    const fLink = (source, target) => {
      let vx = target.x - source.x
      let vy = target.y - source.y
      const vr = Math.sqrt(vx * vx + vy * vy)
      const sin = vy / vr
      const cos = vx / vr
      const sx = vy ? Math.max(-1, Math.min(cos / sin, 1)) * Math.sign(sin) : Math.sign(vx)
      const sy = vx ? Math.max(-1, Math.min(sin / cos, 1)) * Math.sign(cos) : Math.sign(vy)
      const targetR = (f) => Math.abs(target.r) * f
      let offsetxt, offsetyt, offsetxs, offsetys
      if (source.r > 0) {
        // offset from edge of circle
        offsetxs = cos * source.r
        offsetys = sin * source.r
      } else {
        // offset from edge of square
        offsetxs = -sx * source.r * whRatio
        offsetys = -sy * source.r
      }
      if (target.r > 0) {
        // offset to edge of circle
        offsetxt = cos * targetR(1)
        offsetyt = sin * targetR(1)
      } else {
        // offset to edge of square
        offsetxt = sx * targetR(whRatio)
        offsetyt = sy * targetR(1)
      }
      source = {x: source.x + offsetxs, y: source.y + offsetys}
      target = {x: target.x - offsetxt, y: target.y - offsetyt}
      vx = target.x - source.x
      vy = target.y - source.y
      const line = d3.line().x(d => d[0]).y(d => d[1])
        .curve(d3.curveBasis)
      const points = [[source.x, source.y]]
      if (Math.abs(vx) > Math.abs(vy)) {
        target.x = target.x - markerWidth * 2 * Math.sign(vx)
        // points.push([source.x + vx / 2, source.y])
        // points.push([source.x + vx / 2, target.y])
      } else {
        target.y = target.y - markerWidth * 2 * Math.sign(vy)
        // points.push([source.x, source.y + vy / 2])
        // points.push([target.x, source.y + vy / 2])
      }
      points.push([target.x, target.y])
      return line(points)
    }

    const zoomTo = (v) => {
      if (!node) return

      // for elementary tasks, set radius to 1/6
      const k = zoomF(v[2])
      const rk = (node) => {
        return node.id === focus.id && !isWorkflow(node) ? 1 / 6 : 1
      }
      const bound = (d) => {
        let x = (d.x - v[0]) * k + width / 2
        let y = (d.y - v[1]) * k + height / 2
        const s = 0.8 * scaling
        if (d.type === 'tasks' && d.id !== source.id) {
          x = Math.min(Math.max(x, s * width), (1 - s) * width)
          y = Math.min(Math.max(y, s * height), (1 - s) * height)
        }
        return [x, y]
      }

      view = v
      node
        .attr('transform', d => {
          const [x, y] = bound(d)
          return `translate(${x},${y})`
        })
      node.selectAll('circle').attr('r', d => d.r * rk(d) * k)
      node.selectAll('rect')
        .attr('x', d => -d.r * rk(d) * k * whRatio)
        .attr('y', d => -d.r * rk(d) * k)
        .attr('rx', d => d.r * rk(d) * k * 0.05)
        .attr('width', d => d.r * 2 * rk(d) * k * whRatio)
        .attr('height', d => d.r * 2 * rk(d) * k)
      line.attr('d', d => {
        const translate = (node) => {
          const isCircle = ['inputs', 'outputs', 'inputs-outputs'].includes(node.type)
          const dr = node.r * rk(node) * k
          const [x, y] = bound(node)
          return {
            x: x,
            y: y,
            r: isCircle ? dr : -dr
          }
        }
        return fLink(translate(d.source), translate(d.target))
      })
      node.selectAll('text').attr('y', d => -1.2 * d.r * rk(d) * k)
    }

    const zoom = (d) => {
      focus = d
      if (!focus) return
      const transition = svg.transition()
        .duration(500)
        .tween('zoom', d => {
          const i = d3.interpolateZoom(view, [focus.x, focus.y, focus.r * 2])
          return t => zoomTo(i(t))
        })

      node.selectAll('text')
        .transition(transition)
    }

    const update = (source) => {
      if (!source) return

      if (!source.nodes || !source.nodes.length) return

      dragged = false

      const setShowNodes = (inputValue) => {
        const value = []
        if (hasError) setCurrentNode(currentNode)
        let sourceNode = currentNode
        if (['inputs', 'outputs'].includes(currentNode.type)) {
          sourceNode = currentNode.parent
        }
        if (!sourceNode) return
        let total = 1
        const nodes = sourceNode.nodes.filter(node => node.type === (currentNode.type || 'tasks'))
        if (nodes.length) total = nodes[0].total
        const toNumber = (value) => {
          if (value.includes('%')) return Math.floor(parseFloat(value.replace('%')) * total / 100)
          if (value.includes('.')) return Math.floor(parseFloat(value.replace('%')) * total)
          value = parseInt(value)
          if (value < 0) value = total + value
          return value
        }
        for (const query of inputValue.split(',')) {
          try {
            const q = query.split(':')
            if (q.length === 2) {
              value.push(...range(...q.map((qi, n) => qi ? toNumber(qi) : n === 0 ? 0 : total)))
            } else {
              value.push(toNumber(q[0]))
            }
          } catch (error) {
            console.error(error)
            hasError = true
            setCurrentNode(null)
            break
          }
        }
        if (!value.length) return
        dquery['resolveIndices'][currentNode.type || 'tasks'] = [...new Set(value)]
        source.nodes = source.nodes.filter(node => node.id !== currentNode.id)
        // set children and links to recalculate node positions and links
        const d = source
        d.children = undefined
        d.links = undefined
        d.parent = source.parent
        resolveSection(d, query).then(d => {
          if (!d) return
          getNodes(d, dquery).then(nodes => {
            d.nodes = nodes
            getLinks(d, query).then(links => {
              d.links = links
              update(d)
              zoom(d)
            }).catch(asyncError)
          }).catch(asyncError)
          .catch(asyncError)
        })
      }

      d3.select('#nodes-filter-clear')
        .on('click', () => {
          const mid = (query.maxNodes || 6) / 2
          setShowNodes(`:${mid},${-mid}:`)
      })

      d3.select('#nodes-filter-enter')
        .on('keydown', () => {
          const inputValue = d3.event.target.value
          if (d3.event.key === 'Enter' && inputValue) setShowNodes(inputValue)
        })

      d3.select('#backbutton')
        .on('click', () => {
          previousNode = source.parent
          handleClickIcon(source)
        })

      d3.select('#resetbutton')
        .on('click', () => {
          previousNode = 'root'
          root.nodes = null
          root.links = null
          root.children = null
          dquery.resolveIndices = {}
          resolveSection(root, query).then(d => {
            if (!d) return
            getNodes(d, query).then(nodes => {
              d.nodes = nodes
              getLinks(d, query).then(links => {
                d.links = links
                update(d)
                zoom(d)
              }).catch(asyncError)
            }).catch(asyncError)
            .catch(asyncError)
          })
          svg.call(zoomBehaviors.transform, d3.zoomIdentity)
          if (enableForce) simulation.alphaTarget(0)
        })

      d3.select('#legendtogglebutton')
        .on('click', () => {
          const visibility = legend.attr('visibility') === 'visible' ? 'hidden' : 'visible'
          legend.attr('visibility', visibility)
          setShowLegend(visibility === 'visible')
        })

      d3.select('#forcetogglebutton')
        .on('click', () => {
          enableForce = !enableForce
          setEnableForce(enableForce)
          update(source)
          zoom(source)
          if (enableForce) {
            simulation.alphaTarget(0.7).restart()
          } else {
            simulation.alphaTarget(0)
          }
        })

      // set ids and size
      const maxLength = Math.max(...source.nodes
        .filter(node => node.type === 'tasks').map(node => node.nChildren || 0))
      source.size = source.nChildren / maxLength
      source.id = id
      id = id + 1
      source.nodes.forEach((node) => {
        node.size = isWorkflow(node) ? node.nChildren / maxLength : 1
        // put a lower limit on size so node will not get too small
        node.size = Math.max(node.size, 0.1)
        node.id = id
        id = id + 1
      })

      setNodesPosition(source)

      nodes = [source, ...source.nodes]
      const links = source.links

      links.forEach((link) => {
        link.id = `${link.source.id}-${link.target.id}`
      })

      addLinkMarkers(links)

      if (enableForce) {
        const k = zoomF(source.r * 2)
        simulation
          .nodes(source.nodes.filter(node => node.type === 'tasks'))
          .on('tick', () => zoomTo(view))

        simulation
          .force('charge')
          .strength(-10 / k ** 2)

        simulation
          .force('link')
          .strength(0.005 / k)
          .distance(source.r)
          .links(source.links.filter(link => link.source.type === 'tasks' && link.target.type === 'tasks'))

        simulation
          .force('center', d3.forceCenter(source.x, source.y))

        simulation
          .force('collide', d3.forceCollide(source.r / 3).strength(2))
      }

      const handleMouseOverIcon = (d) => {
        d3.select(`#icon-${d.id}`).style('stroke-opacity', 1)
        if (d.id === source.id) {
          if (!previousNode || previousNode === 'root') return
          setTooltipContent(<p>Click to go back up</p>)
          return
        }
        if (['inputs', 'outputs'].includes(d.type)) {
          setTooltipContent(<p>Click to switch {d.type} filter</p>)
        }
        const sectionType = d.sectionType === 'tasks' ? 'task' : 'workflow'
        if (d.type === 'tasks') setTooltipContent(<p>Click to expand {sectionType}</p>)
      }

      const handleMouseOutIcon = (d) => {
        setTooltipContent('')
        d3.select(`#icon-${d.id}`).style('stroke-opacity', 0.5)
      }

      const handleClickIcon = (d) => {
        if (dragged) d.children = null
        d3.event.stopPropagation()
        setCurrentNode(d)
        dquery['resolveIndices'][d.type || 'tasks'] = null
        currentNode = d
        if (!d.nodes || !d.nodes.length) {
          return
        }
        if (d.id === source.id) {
          if (!previousNode || previousNode === 'root') return
          setCurrentNode(previousNode)
          currentNode = previousNode
          update(previousNode)
          zoom(previousNode)
          previousNode = previousNode.parent
          return
        }
        d.parent = source
        previousNode = d.parent
        resolveSection(d, query).then(d => {
          if (!d) return
          getNodes(d, query).then(nodes => {
            d.nodes = nodes
            getLinks(d, query).then(links => {
              d.links = links
              update(d)
              zoom(d)
            }).catch(asyncError)
          }).catch(asyncError)
          .catch(asyncError)
        })
      }

      node = gNode.selectAll('g')
        .attr('visibility', 'hidden')
        .attr('pointer-events', 'none')
        .exit().remove()
        .data(nodes, d => d.id)
        .attr('cursor', d => d.nodes ? null : 'pointer')
        .join('g')
        .call(dragBehaviors)

      node
        .filter(d => d.type === 'tasks' || d.id === source.id)
        .append('rect')
        .attr('class', 'icon')
        .attr('id', d => `icon-${d.id}`)
        .attr('stroke', color.outline)
        .attr('fill', d => nodeColor(d))
        .attr('fill-opacity', d => {
          if (d.type === 'link') return 0
          if (!d.nodes || !d.nodes.filter(node => node.type === 'tasks').length) return 1
          if (d.id === source.id) return 0.2
          return 1
        })
        .on('mouseover', handleMouseOverIcon)
        .on('mouseout', handleMouseOutIcon)
        .on('click', handleClickIcon)

      node
        .filter(d => d.type === 'inputs' || d.type === 'outputs' || d.type === 'inputs-outputs')
        .append('circle')
        .attr('class', 'icon')
        .attr('id', d => `icon-${d.id}`)
        .attr('stroke', color.outline)
        .attr('fill', d => nodeColor(d))
        .on('mouseover', handleMouseOverIcon)
        .on('mouseout', handleMouseOutIcon)
        .on('click', handleClickIcon)

      node.append('text')
        .attr('class', 'text')
        .attr('fill', color.text)
        .attr('font-weight', d => d.id === source.id ? 'bold' : 'none')
        .attr('id', d => `text-${d.id}`)
        .text(d => trimName(d.name))
        .style('font-size', d => d.id === nodes[0].id ? 18 : 14)
        .on('click', d => {
          d3.event.stopPropagation()
          if (!d.entryId || !d.parent) return
          let path = `entry/id/${d.entryId}`
          const sectionPath = d.path ? d.path.replace(/\/(?=\d)/g, ':') : null
          path = isWorkflow(d) ? path : sectionPath ? `${path}/data${sectionPath}` : path
          const url = getUrl(path)
          history.push(url)
        })
        .on('mouseover', d => {
          if (!d.type || !d.parent) return
          if (!d.sectionType) return
          d3.select(`#text-${d.id}`).style('font-weight', 'bold')
            .text(d.name)
          const text = isWorkflow(d) ? 'overview page' : 'archive section'
          if (d.entryId) {
            setTooltipContent(<p>Click to go to {text} for entry<br/>{d.entryId}</p>)
          }
        })
        .on('mouseout', d => {
          setTooltipContent('')
          d3.select(`#text-${d.id}`).style('font-weight', null)
            .text(d => trimName(d.name))
        })

      const link = gLink.selectAll('path')
        .attr('visibility', 'hidden')
        .attr('pointer-events', 'none')
        .exit().remove()
        .data(links)

      line = link.enter()
        .append('path')
        .attr('pointer-events', 'all')
        .attr('id', d => `link-${d.id}`)
        .attr('marker-end', d => `url(#${d.id})`)
        .attr('visibility', 'visible')
        .on('click', d => {
          d3.event.stopPropagation()
          const parent = d.source.parent || d.target.parent
          if (!parent) return
          // use the parent node to containt source and target nodes
          const linkNode = parent
          const sourceNode = d.source
          const targetNode = d.target
          if (sourceNode.type !== 'tasks' || targetNode.type !== 'tasks') return

          const store = (node, recursive) => {
            if (typeof node !== 'object') return
            node._name = node.name
            node._type = node.type
            node.children = null
            if (node.links) {
              node._links = [...node.links]
              node.links = null
            }
            if (node.nodes) {
              node._nodes = [...node.nodes]
              if (recursive) {
                node.nodes.forEach(n => {
                  store(n)
                })
              }
            }
          }

          const restore = (node) => {
            if (typeof node !== 'object' || !node._nodes) return
            node.nodes = node._nodes
            node.links = node._links
            node.children = null
            node.name = node._name || node.name
            node.type = node._type || node.type
          }

          // store the original nodes info
          store(linkNode)
          store(sourceNode)
          store(targetNode)

          // include into sourceNode the tagetNode inputs
          const inputs = [...targetNode.nodes.filter(node => node.type === 'inputs')]
          let url = {}
          sourceNode.nodes.forEach(node => {
            url[node.url] = node.name
          })
          inputs.forEach(input => {
            if (Object.keys(url).includes(input.url)) {
              sourceNode.nodes.push(input)
              input._name = input.name
              const sourceName = url[input.url]
              if (sourceName !== input.name) input.name = `${input.name}/${url[input.url]}`
              input._type = input.type
              input.type = 'inputs-outputs'
            }
          })

          linkNode.type = 'link'
          // parent nodes should contain the source, target
          linkNode.nodes = [sourceNode, targetNode]
          // and the targetNode inputs
          linkNode.nodes.push(...inputs)

          // and the sourceNode outputs
          url = inputs.map(node => node.url)
          const outputs = sourceNode.nodes.filter(node => node.type === 'outputs' && !url.includes(node.url))
          linkNode.nodes.push(...outputs)

          // and the tagetNode outputs
          linkNode.nodes.push(...targetNode.nodes.filter(node => node.type === 'outputs'))

          // generate links for this configuration
          getLinks(linkNode, query).then(links => {
            for (const [index, link] of links.entries()) {
              if (link.target.url === sourceNode.url) {
                // flip source and target since node is output of sourceNode
                links[index] = {source: link.target, target: link.source, id: link.id, label: 'Output'}
              } else if (link.target.url === linkNode.url || link.source.url === linkNode.url) {
                // remove link to and from rootNode
                links[index] = null
              } else if (link.source.url === sourceNode.url && link.target.url === targetNode.url) {
                // remove link from sourceNode to targetNode, unnecessary
                links[index] = null
              }
            }
            linkNode.links = links.filter(link => link)
            linkNode.name = null
            // linkNode.type = 'link'
            update(linkNode)
            zoom(linkNode)
            previousNode = linkNode
            // reset original nodes info
            inputs.forEach(input => {
              input.type = 'inputs'
              input.name = input._name || input.name
              input.children = null
            })
            linkNode.type = linkNode._type
            restore(linkNode)
            restore(sourceNode)
            restore(targetNode)
            setNodesPosition(linkNode)
        })
        })
        .on('mouseover', d => {
          d3.select(`#link-${d.id}`).style('stroke-opacity', 1.0)
          svg.select(`.marker-${d.id}`).attr('fill-opacity', 1.0)
          d3.select(`#icon-${d.source.id}`).style('stroke', color.link).style('stroke-opacity', 1.0)
          d3.select(`#icon-${d.target.id}`).style('stroke', color.link).style('stroke-opacity', 1.0)
          setTooltipContent(<p>{d.label}</p>)
        })
        .on('mouseout', d => {
          d3.select(`#link-${d.id}`).style('stroke-opacity', 0.5)
          svg.select(`.marker-${d.id}`).attr('fill-opacity', 0.5)
          d3.select(`#icon-${d.source.id}`).style('stroke', color.outline).style('stroke-opacity', 0.5)
          d3.select(`#icon-${d.target.id}`).style('stroke', color.outline).style('stroke-opacity', 0.5)
          setTooltipContent('')
        })
    }

    resolveSection(source, query).then(source => {
      if (!source) return
      getNodes(source, query).then(nodes => {
        source.nodes = nodes
        source.x = 0
        source.y = 0
        source.r = width / 2
        getLinks(source, query).then(links => {
          source.links = links
          focus = source
          root = source
          setCurrentNode(source)
          currentNode = source
          update(source)
          zoomTo([source.x, source.y, source.r * 2])
        }).catch(asyncError)
      }).catch(asyncError)
      .catch(asyncError)
    })
  }, [
    history,
    setTooltipContent,
    setCurrentNode,
    setShowLegend,
    setEnableForce,
    query,
    source,
    svgRef,
    finalLayout,
    asyncError
  ])

  return <svg className={classes.root} ref={svgRef}></svg>
})

Graph.propTypes = {
  source: PropTypes.object.isRequired,
  query: PropTypes.object.isRequired,
  layout: PropTypes.object,
  setTooltipContent: PropTypes.any,
  setCurrentNode: PropTypes.any,
  setShowLegend: PropTypes.any,
  setEnableForce: PropTypes.any
}

const WorkflowCard = React.memo(({archive}) => {
  const {api} = useApi()
  const [tooltipContent, setTooltipContent] = useState('')
  const [tooltipPosition, setTooltipPosition] = useState({x: undefined, y: undefined})
  const [showLegend, setShowLegend] = useState(true)
  const [enableForce, setEnableForce] = useState(false)
  const [currentNode, setCurrentNode] = useState({type: 'tasks'})
  const [inputValue, setInputValue] = useState('')
  const query = useMemo(() => ({
    api: api,
    required: { 'workflow2': '*', 'metadata': '*', 'data': '*' },
    sectionKeys: ['inputs', 'tasks', 'outputs'],
    maxNodes: 6
  }), [api])

  const graph = useMemo(() => {
    if (!archive || !archive.workflow2) return ''

    const source = {
      path: '/workflow2',
      archive: archive
    }

    return <Graph
      source={source}
      query={query}
      setTooltipContent={setTooltipContent}
      setCurrentNode={setCurrentNode}
      setShowLegend={setShowLegend}
      setEnableForce={setEnableForce}
    ></Graph>
  }, [archive, query])

  let nodeType = 'tasks'
  let nodesCount = 0
  if (currentNode) {
    nodeType = currentNode.type || 'tasks'
    let sourceNode = currentNode
    if (['inputs', 'outputs'].includes(nodeType)) sourceNode = currentNode.parent || currentNode
    const nodes = sourceNode.nodes ? sourceNode.nodes.filter(node => node.type === nodeType) : []
    if (nodes.length) nodesCount = nodes[0].total
  }

  let label = `No ${nodeType.slice(0, -1)} to show`
  if (nodesCount) label = `Filter ${nodeType} (N=${nodesCount})`

  const actions = <div>
    <FormControl margin='none'>
      <TextField
        error={!currentNode}
        size='medium'
        id='nodes-filter-enter'
        label={currentNode ? label : 'Invalid input'}
        disabled={!nodesCount || !currentNode}
        placeholder={'Enter range: 5, 4:-1, :50%'}
        variant='outlined'
        value={inputValue}
        onChange={(event) => setInputValue(event.target.value)}
        InputProps={{
          endAdornment: (
            <IconButton id='nodes-filter-clear' onClick={() => setInputValue('')} size='medium'>
              <Clear />
            </IconButton>
          )
        }}
      ></TextField>
    </FormControl>
    <IconButton id='forcetogglebutton'>
      <Tooltip title={`${enableForce ? 'Disable' : 'Enable'} force simulation`}>
        {enableForce ? <StopSharp /> : <PlayArrowSharp />}
      </Tooltip>
    </IconButton>
    <IconButton id='legendtogglebutton'>
      <Tooltip title={showLegend ? 'Hide legend' : 'Show legend'}>
        {showLegend ? <LabelOff /> : <Label />}
      </Tooltip>
    </IconButton>
    <IconButton id='backbutton'>
      <Tooltip title="Back">
        <Undo />
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
