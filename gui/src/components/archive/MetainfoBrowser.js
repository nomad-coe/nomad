/*
 * Copyright The NOMAD Authors.
 *
 * This file is part of NOMAD. See https://nomad-lab.eu for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
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
import React, { useMemo, useEffect, useRef, useLayoutEffect, useContext, useState } from 'react'
import PropTypes from 'prop-types'
import { useRecoilValue, useRecoilState, atom } from 'recoil'
import { configState } from './ArchiveBrowser'
import Browser, { Item, Content, Compartment, Adaptor, laneContext, formatSubSectionName } from './Browser'
import { Typography, Box, makeStyles, Grid, FormGroup, TextField, Button, Tooltip, Link } from '@material-ui/core'
import { metainfoDef, resolveRef, vicinityGraph, rootSections, path as metainfoPath, packagePrefixes, defsByName } from './metainfo'
import * as d3 from 'd3'
import blue from '@material-ui/core/colors/blue'
import teal from '@material-ui/core/colors/teal'
import lime from '@material-ui/core/colors/lime'
import purple from '@material-ui/core/colors/purple'
import grey from '@material-ui/core/colors/grey'
import Markdown from '../Markdown'
import { JsonCodeDialogButton } from '../buttons/CodeDialogButton'
import Histogram from '../Histogram'
import { appBase } from '../../config'
import { useHistory, useRouteMatch } from 'react-router-dom'
import Autocomplete from '@material-ui/lab/Autocomplete'
import { useApi } from '../api'
import { useErrors } from '../errors'

export const help = `
The NOMAD *metainfo* defines all quantities used to represent archive data in
NOMAD. You could say it is the archive *schema*. You can browse this schema and
all its definitions here.

The NOMAD metainfo contains three different *kinds* of definitions:

- **sections**: A section is a nested groups of quantities that allow a hierarchical data structure
- **values**: Actual quantities that contain data
- **references**: References that allow to connect related sections.

All definitions have a name that you can search for. Furthermore, all definitions
are organized in packages. There is a *common* pkg with definitions that are
used by all codes and there are packages for each code with code specific definitions.
You can select the pkg to browse below.

Depending on the selected pkg, there are quite a large number of definitions.
You can use the *definition* field to search based on definition names.

All definitions are represented as *cards* below. Click on the various card items
to expand sub-sections, open values or references, hide and show compartments, or
collapse cards again. The highlighted *main* card cannot be collapsed. The
shapes in the background represent section containment (grey) and
reference (blue) relations.

If you bookmark this page, you can save the definition represented by the highlighted
*main* card.

To learn more about the meta-info, visit the [meta-info documentation](${appBase}/docs/metainfo.html).
`

const showInnerSectionDefinitions = false

function defCompare(a, b) {
  return a.name.localeCompare(b.name)
}

export const metainfoConfigState = atom({
  key: 'metainfoConfig',
  default: {
    'packagePrefix': 'nomad'
  }
})

export function MetainfoPage() {
  return <Box margin={3}>
    <MetainfoBrowser />
  </Box>
}

export default function MetainfoBrowser() {
  return <Browser
    adaptor={new MetainfoRootAdaptor()}
    form={<MetainfoConfigForm />}
  />
}

function MetainfoConfigForm(props) {
  const [config, setConfig] = useRecoilState(metainfoConfigState)

  const history = useHistory()
  const { url } = useRouteMatch()

  const searchOptions = useMemo(() => {
    return Object.keys(defsByName).reduce((results, name) => {
      const defsForName = defsByName[name].filter(def => (
        !def.extends_base_section && def.m_def !== 'SubSection' && (
          def._package.name.startsWith(config.packagePrefix) || def._package.name.startsWith('nomad'))
      ))
      results.push(...defsForName.map(def => {
        let label = def.name
        if (defsForName.length > 1) {
          if (def._parentSections?.length) {
            label = `${label} (${def._parentSections[0].name})`
          } else {
            label = `${label} (${def._qualifiedName.split(':')[0]})`
          }
        }
        return {
          path: url + '/' + metainfoPath(def),
          label: `${label} [${def.m_def.toLowerCase()}]`
        }
      }))
      return results
    }, []).sort((a, b) => a.label.localeCompare(b.label))
  }, [config.packagePrefix, url])

  return (
    <Box marginTop={-2}>
      <FormGroup row style={{alignItems: 'flex-end'}}>
        <Autocomplete
          options={searchOptions}
          getOptionLabel={(option) => option.label}
          style={{ width: 500 }}
          onChange={(_, value) => {
            if (value) {
              history.push(value.path)
            }
          }}
          renderInput={(params) => <TextField {...params} label="search" margin="normal" />}
        />
        <Box margin={1} />
        <Autocomplete
          value={config.packagePrefix}
          options={Object.keys(packagePrefixes)}
          getOptionLabel={option => option.replace(/parser/g, '')}
          style={{ width: 350 }}
          onChange={(_, value) => setConfig({...config, packagePrefix: value})}
          renderInput={(params) => <TextField {...params} label="source" margin="normal" />}
        />
      </FormGroup>
    </Box>
  )
}

export function metainfoAdaptorFactory(obj) {
  if (obj.m_def === 'Section') {
    return new SectionDefAdaptor(obj)
  } else if (obj.m_def === 'SubSection') {
    return new SubSectionDefAdaptor(obj)
  } else if (obj.m_def === 'Quantity') {
    return new QuantityDefAdaptor(obj)
  } else if (obj.m_def === 'Category') {
    return new CategoryDefAdaptor(obj)
  } else {
    throw new Error('Unknown metainfo definition type')
  }
}

class MetainfoAdaptor extends Adaptor {
  itemAdaptor(key) {
    if (key === '_reference') {
      return metainfoAdaptorFactory(resolveRef(this.e.type.type_data))
    } else if (key.startsWith('_category:')) {
      const categoryName = key.split(':')[1]
      return metainfoAdaptorFactory(this.e.categories.map(ref => resolveRef(ref)).find(categoryDef => categoryDef.name === categoryName))
    } else if (key === '_metainfo') {
      return metainfoAdaptorFactory(metainfoDef(this.e.m_def))
    } else if (this.e[key]) {
      return metainfoAdaptorFactory(resolveRef(this.e[key]))
    } else {
      throw new Error('Unknown item key')
    }
  }
}

export class MetainfoRootAdaptor extends MetainfoAdaptor {
  itemAdaptor(key) {
    const rootSection = rootSections.find(def => def.name === key)
    if (rootSection) {
      return metainfoAdaptorFactory(rootSection)
    } else if (packagePrefixes[key]) {
      return new PackagePrefixAdaptor(packagePrefixes[key])
    } else {
      super.itemAdaptor(key)
    }
  }
  render() {
    return <Metainfo />
  }
}

export class PackagePrefixAdaptor extends MetainfoAdaptor {
  itemAdaptor(key) {
    const [type, name] = key.split('@')
    const def = Object.keys(this.e)
      .map(key => this.e[key])
      .reduce((value, pkg) => value || pkg[type].find(def => def._qualifiedName === name), null)
    return metainfoAdaptorFactory(def)
  }
  render() {
    const sectionDefs = Object.keys(this.e)
      .map(key => this.e[key])
      .reduce((defs, pkg) => {
        pkg.section_definitions.forEach(def => defs.push(def))
        return defs
      }, [])
    const categoryDefs = Object.keys(this.e)
      .map(key => this.e[key])
      .reduce((defs, pkg) => {
        pkg.category_definitions.forEach(def => defs.push(def))
        return defs
      }, [])

    return <Content>
      <Compartment title="Sections">
        {sectionDefs.filter(def => !def.extends_base_section).sort(defCompare).map(def => {
          const key = `section_definitions@${def._qualifiedName}`
          return <Item key={key} itemKey={key}>
            <Typography>{def.name}</Typography>
          </Item>
        })}
      </Compartment>
      <Compartment title="Section Extensions">
        {sectionDefs.filter(def => def.extends_base_section).sort(defCompare).map(def => {
          const key = `section_definitions@${def._qualifiedName}`
          return <Item key={key} itemKey={key}>
            <Typography>{def.name}</Typography>
          </Item>
        })}
      </Compartment>
      <Compartment title="Categories">
        {categoryDefs.sort(defCompare).map(def => {
          const key = `category_definitions@${def._qualifiedName}`
          return <Item key={key} itemKey={key}>
            <Typography>{def.name}</Typography>
          </Item>
        })}
      </Compartment>
    </Content>
  }
}

function Metainfo(props) {
  return <Content>
    <Compartment title="archive root section">
      <Item itemKey="EntryArchive">
        <Typography>EntryArchive</Typography>
      </Item>
    </Compartment>
    <Compartment title="other root sections">
      {rootSections.filter(def => def.name !== 'EntryArchive').map(def => (
        <Item key={def.name} itemKey={def.name}>
          <Typography>
            {def.name}
          </Typography>
        </Item>
      ))}
    </Compartment>
    <Compartment title="sources">
      {Object.keys(packagePrefixes).map(key => <Item key={key} itemKey={key}>
        <Typography>{key.replace(/parser$/, '')}</Typography>
      </Item>)}
    </Compartment>
  </Content>
}

export class SectionDefAdaptor extends MetainfoAdaptor {
  itemAdaptor(key) {
    if (key === '_baseSection') {
      return metainfoAdaptorFactory(resolveRef(this.e.base_sections[0]))
    }

    if (key.includes('@')) {
      const [type, name] = key.split('@')
      if (type === 'innerSectionDef') {
        const innerSectionDef = this.e.inner_section_definitions.find(def => def.name === name)
        if (innerSectionDef) {
          return metainfoAdaptorFactory(innerSectionDef)
        }
      }
    }

    const property = this.e._properties[key]
    if (property) {
      return metainfoAdaptorFactory(property)
    }

    return super.itemAdaptor(key)
  }
  render() {
    return <SectionDef def={this.e} />
  }
}

class SubSectionDefAdaptor extends MetainfoAdaptor {
  constructor(e) {
    super(e)
    this.sectionDefAdaptor = new SectionDefAdaptor(resolveRef(e.sub_section))
  }
  itemAdaptor(key) {
    return this.sectionDefAdaptor.itemAdaptor(key)
  }
  render() {
    return <SubSectionDef def={this.e} />
  }
}

class QuantityDefAdaptor extends MetainfoAdaptor {
  render() {
    return <QuantityDef def={this.e} />
  }
}

class CategoryDefAdaptor extends MetainfoAdaptor {
  render() {
    return <Content>
      <Definition def={this.e} />
      <DefinitionDetails def={this.e} />
    </Content>
  }
}

function SectionDefContent({def}) {
  const config = useRecoilValue(configState)
  const metainfoConfig = useRecoilValue(metainfoConfigState)
  const filter = def.extends_base_section ? () => true : def => {
    if (def._package.name.startsWith('nomad')) {
      return true
    }
    if (def._package.name.startsWith('nexus')) {
      return true
    }
    if (metainfoConfig.packagePrefix) {
      return def._package.name.startsWith(metainfoConfig.packagePrefix)
    }
    if (config.showCodeSpecific) {
      return true
    }
    return false
  }

  return <React.Fragment>
    <DefinitionProperties def={def} />
    {def.base_sections.length > 0 &&
      <Compartment title="base section">
        {def.base_sections.map(baseSectionRef => {
          const baseSection = resolveRef(baseSectionRef)
          return <Item key={baseSectionRef} itemKey="_baseSection">
            <Typography>{baseSection.name}</Typography>
          </Item>
        })}
      </Compartment>
    }
    <Compartment title="sub section definitions">
      {def.sub_sections.filter(filter)
        .map(subSectionDef => {
          const key = subSectionDef.name
          const categories = subSectionDef.categories?.map(c => resolveRef(c))
          const unused = categories?.find(c => c.name === 'Unused')
          return <Item key={key} itemKey={key}>
            <Typography component="span" color={unused && 'error'}>
              <Box fontWeight="bold" component="span">
                {formatSubSectionName(subSectionDef.name)}
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
          const categories = quantityDef.categories?.map(c => resolveRef(c))
          const unused = categories?.find(c => c.name === 'Unused')
          return <Item key={key} itemKey={key}>
            <Box component="span" whiteSpace="nowrap">
              <Typography component="span" color={unused && 'error'}>
                <Box fontWeight="bold" component="span">
                  {quantityDef.name}
                </Box>
              </Typography>
            </Box>
          </Item>
        })
      }
    </Compartment>
    {showInnerSectionDefinitions && <Compartment title="inner section definitions">
      {def.inner_section_definitions.filter(filter)
        .map(innerSectionDef => {
          const key = `innerSectionDef@${innerSectionDef.name}`
          const categories = innerSectionDef.categories?.map(c => resolveRef(c))
          const unused = categories?.find(c => c.name === 'Unused')
          return <Item key={key} itemKey={key}>
            <Box component="span" whiteSpace="nowrap">
              <Typography component="span" color={unused && 'error'}>
                <Box fontWeight="bold" component="span">
                  {innerSectionDef.name}
                </Box>
              </Typography>
            </Box>
          </Item>
        })
      }
    </Compartment>}
    <DefinitionDetails def={def} />
  </React.Fragment>
}
SectionDefContent.propTypes = ({
  def: PropTypes.object
})

function SectionDef({def}) {
  return <Content style={{backgroundColor: 'grey'}}>
    <Definition def={def} kindLabel="section definition" />
    <SectionDefContent def={def} />
  </Content>
}
SectionDef.propTypes = ({
  def: PropTypes.object
})

function SubSectionDef({def}) {
  const sectionDef = resolveRef(def.sub_section)
  return <React.Fragment>
    <Content>
      <Title def={def} isDefinition kindLabel="sub section definition" />
      <DefinitionDocs def={sectionDef} />
      <SectionDefContent def={sectionDef} />
    </Content>
  </React.Fragment>
}
SubSectionDef.propTypes = ({
  def: PropTypes.object
})

function DefinitionProperties({def, children}) {
  if (!(children || def.aliases?.length || def.deprecated || Object.keys(def.more).length)) {
    return ''
  }

  return <Compartment title="properties">
    {children}
    {def.aliases?.length && <Typography><b>aliases</b>:&nbsp;{def.aliases.map(a => `"${a}"`).join(', ')}</Typography>}
    {def.deprecated && <Typography><b>deprecated</b>: {def.deprecated}</Typography>}
    {Object.keys(def.more).map((moreKey, i) => (
      <Typography key={i}><b>{moreKey}</b>:&nbsp;{String(def.more[moreKey])}</Typography>
    ))}
  </Compartment>
}
DefinitionProperties.propTypes = ({
  def: PropTypes.object,
  children: PropTypes.any
})

function QuantityDef({def}) {
  return <Content>
    <Definition def={def} kindLabel="quantity definition"/>
    <DefinitionProperties def={def}>
      {def.type.type_kind !== 'reference'
        ? <Typography>
          <b>type</b>:&nbsp;
          {Array.isArray(def.type.type_data) ? def.type.type_data.join(', ') : def.type.type_data}&nbsp;
          {def.type.type_kind !== 'data' && `(${def.type.type_kind})`}
        </Typography>
        : <Item itemKey="_reference">
          <Typography>
            <b>referenced section</b>
          </Typography>
        </Item>}
      <Typography>
        <b>shape</b>:&nbsp;
        [{def.shape.join(', ')}]
      </Typography>
      {def.unit &&
        <Typography><b>unit</b>:&nbsp;{def.unit}</Typography>}
      {def.derived && <Typography><b>derived</b></Typography>}
    </DefinitionProperties>
  </Content>
}
QuantityDef.propTypes = ({
  def: PropTypes.object
})

function DefinitionDocs({def}) {
  return <React.Fragment>
    {def.description && !def.extends_base_section &&
      <Compartment title="description">
        <Box marginTop={1} marginBottom={1}>
          {def._qualifiedName.startsWith('nexus')
            ? <Typography>{def.description}</Typography>
            : <Markdown>{def.description}</Markdown>}
        </Box>
      </Compartment>
    }
    {def.links?.length &&
      <Compartment title="links">
        {def.links.map((link, i) => <div key={i}>
          <Link href={link} key={i}>{link.includes('manual.nexusformat.org') ? 'nexus manual' : link}</Link>
        </div>)}
      </Compartment>
    }
  </React.Fragment>
}
DefinitionDocs.propTypes = {
  def: PropTypes.object.isRequired
}

function Definition({def, ...props}) {
  return <React.Fragment>
    <Title def={def} isDefinition {...props} />
    <DefinitionDocs def={def} />
  </React.Fragment>
}
Definition.propTypes = {
  def: PropTypes.object.isRequired
}

function DefinitionDetails({def, ...props}) {
  const {api} = useApi()
  const {raiseError} = useErrors()
  const lane = useContext(laneContext)
  const isLast = !lane.next
  const [usage, setUsage] = useState(null)
  const [showUsage, setShowUsage] = useState(false)

  const quantityPath = useMemo(() => {
    const path = lane.path.split('/')
    const index = path.indexOf('EntryArchive')
    if (index >= 0) {
      return path.slice(index + 1).join('.')
    }
    return null
  }, [lane])

  useEffect(() => {
    if (showUsage) {
      api.post('/entries/query', {
        owner: 'visible',
        query: {
          'quantities:any': [quantityPath]
        },
        aggregations: {
          program_names: {
            terms: {
              quantity: 'results.method.simulation.program_name',
              pagination: {
                page_size: 100 // make sure we get all codes
              }
            }
          }
        }
      }).then(response => {
        const aggData = response.aggregations.program_names.terms.data
        setUsage(aggData)
      }).catch(raiseError)
    }
  }, [api, raiseError, showUsage, quantityPath, setUsage])

  return <React.Fragment>
    {def.categories && def.categories.length > 0 && <Compartment title="Categories">
      {def.categories.map(categoryRef => {
        const categoryDef = resolveRef(categoryRef)
        return <Item key={categoryRef} itemKey={'_category:' + categoryDef.name}>
          <Typography>{categoryDef.name}</Typography>
        </Item>
      })}
    </Compartment>}
    {isLast && !def.extends_base_section && def.name !== 'EntryArchive' &&
      <Compartment title="graph">
        <VicinityGraph def={def} />
      </Compartment>
    }
    {quantityPath &&
      <Compartment title="usage">
        {!showUsage && <Button fullWidth variant="outlined" onClick={() => setShowUsage(true)}>Show usage</Button>}
        {showUsage && !usage && <Typography><i>loading ...</i></Typography>}
        {usage && usage.length > 0 && (
          <Histogram
            data={usage.map(use => ({
              key: use.value,
              name: use.value,
              value: use.count
            }))}
            initialScale={0.5}
            title="Metadata use per code"
          />
        )}
        {usage && usage.length === 0 && (
          <Typography color="error"><i>This metadata is not used at all.</i></Typography>
        )}
      </Compartment>
    }
  </React.Fragment>
}
DefinitionDetails.propTypes = {
  def: PropTypes.object.isRequired
}

const definitionLabels = {
  'Section': 'section',
  'Quantity': 'quantity',
  'SubSection': 'sub section',
  'Category': 'category'
}
export function Title({def, isDefinition, data, kindLabel}) {
  const color = isDefinition ? 'primary' : 'initial'
  return <Compartment>
    <Grid container justifyContent="space-between" wrap="nowrap" spacing={1}>
      <Grid item>
        <Tooltip title={def._qualifiedName || def.name}>
          <Typography color={color} variant="h6">{def.name}</Typography>
        </Tooltip>
        <DefinitionLabel def={def} isDefinition={isDefinition} variant="caption" color={color} />
      </Grid>
      <Grid item>
        <JsonCodeDialogButton
          tooltip={`Show ${(kindLabel + ' ') || ' '}data as JSON`}
          title={`Underlying ${(kindLabel + ' ') || ' '}data as JSON`}
          buttonProps={{size: 'small'}} json={data || def}
        />
      </Grid>
    </Grid>
  </Compartment>
}
Title.propTypes = ({
  def: PropTypes.object.isRequired,
  data: PropTypes.any,
  isDefinition: PropTypes.bool,
  kindLabel: PropTypes.string
})

export function DefinitionLabel({def, isDefinition, ...props}) {
  let label = definitionLabels[def.m_def]
  if (def.extends_base_section) {
    label += ' extension'
  }
  return <Typography {...props}>
    {label}{isDefinition ? ' definition' : ''}
  </Typography>
}
DefinitionLabel.propTypes = ({
  def: PropTypes.object.isRequired,
  isDefinition: PropTypes.bool
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
    const px = w < 400 ? (400 - w) / 2 : 0
    x1 -= px; x2 += px
    graph.viewBox = `${x1} ${y1} ${x2 - x1} ${y2 - y1}`
    graph.aspectRatio = Math.abs((x2 - x1) / (y2 - y1))

    return graph
  }, [def])

  const classes = useVicinityGraphStyles()
  const svgRef = useRef()
  const history = useHistory()

  useEffect(() => {
    const svg = d3.select(svgRef.current)

    const link = svg.select('.links')
      .selectAll('line')
      .data(graph.links)
      .enter().append('line')
      .attr('marker-end', d => `url(#arrowhead-${d.def.m_def || '_none'})`)
      .attr('stroke', d => linkColors[d.def.m_def || '_none'])
      .on('click', d => {
        if (d.def.m_def === 'Quantity') {
          const path = metainfoPath(d.def)
          if (path) {
            history.push(`/metainfo/${path}`)
          }
        }
      })

    const node = svg.select('.nodes')
      .selectAll('g')
      .data(graph.nodes)
      .enter().append('g')

    node.append('circle')
      .attr('r', 10)
      .attr('fill', d => nodeColors[d.def.m_def] || '#000')
      .on('click', d => {
        const path = metainfoPath(d.def)
        if (path) {
          history.push(`/metainfo/${path}`)
        }
      })
      .call(d3.drag()
        .on('drag', d => {
          d.x = d3.event.x
          d.y = d3.event.y
          ticked()
        }))

    node.append('text')
      .text(d => d.def.name)
      .attr('text-anchor', 'middle')
      .attr('y', d => ((d.i % 2) === 0) ? 20 : -14)

    function ticked() {
      link
        .attr('x1', d => d.source.x)
        .attr('y1', d => d.source.y)
        .attr('x2', d => d.target.x)
        .attr('y2', d => d.target.y)

      node
        .attr('transform', d => 'translate(' + d.x + ',' + d.y + ')')
    }
    ticked()
  }, [graph, history, linkColors, nodeColors, svgRef])

  useLayoutEffect(() => {
    svgRef.current.style.height = svgRef.current.clientWidth / graph.aspectRatio
  }, [graph, svgRef])

  return <svg
    className={classes.root} ref={svgRef}
    viewBox={graph.viewBox}
  >
    <g className='links' />
    <g className='nodes' />
    {Object.keys(linkColors).map(colorKey => {
      const color = linkColors[colorKey]
      return <marker
        key={colorKey}
        id={`arrowhead-${colorKey}`}
        viewBox='-0 -5 10 10'
        refX='21'
        refY='0'
        orient='auto'
        markerWidth='3'
        markerHeight='4'
        xoverflow='visible'
      >
        <path d='M 0,-5 L 10 ,0 L 0,5' fill={color} stroke={color} />
      </marker>
    })}
  </svg>
}
VicinityGraph.propTypes = ({
  def: PropTypes.object.isRequired
})
