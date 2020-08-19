
import React, { useMemo, useContext } from 'react'
import PropTypes from 'prop-types'
import { atom, useRecoilState, useRecoilValue } from 'recoil'
import { Box, FormGroup, FormControlLabel, Checkbox, TextField, Typography, makeStyles } from '@material-ui/core'
import { useRouteMatch, useHistory } from 'react-router-dom'
import Autocomplete from '@material-ui/lab/Autocomplete'
import Browser, { Item, Content, Compartment, List, Adaptor } from './Browser'
import { resolveRef, rootSections } from './metainfo'
import { Title, metainfoAdaptorFactory, DefinitionLabel } from './MetainfoBrowser'
import { Matrix, Number } from './visualizations'
import Structure from '../visualization/Structure'
import { StructureViewer } from '@lauri-codes/materia'
import Markdown from '../Markdown'

export const configState = atom({
  key: 'config',
  default: {
    'showMeta': false,
    'showCodeSpecific': false,
    'showAllDefined': false
  }
})

// Context for sharing data that is expensive to create
const ArchiveContext = React.createContext()
const viewer = new StructureViewer()

export default function ArchiveBrowser({data}) {
  const searchOptions = useMemo(() => archiveSearchOptions(data), [data])
  return (
    <ArchiveContext.Provider value={{
      structure: viewer
    }}>
      <Browser
        adaptor={archiveAdaptorFactory(data)}
        form={<ArchiveConfigForm searchOptions={searchOptions} />}
      />
    </ArchiveContext.Provider>
  )
}
ArchiveBrowser.propTypes = ({
  data: PropTypes.object.isRequired
})

function ArchiveConfigForm({searchOptions}) {
  const [config, setConfig] = useRecoilState(configState)
  const handleConfigChange = event => {
    const changes = {[event.target.name]: event.target.checked}
    if (changes.showCodeSpecific) {
      changes.showAllDefined = !changes.showCodeSpecific
    } else if (changes.showAllDefined) {
      changes.showCodeSpecific = !changes.showAllDefined
    }
    setConfig({...config, ...changes})
  }

  const history = useHistory()
  const { url } = useRouteMatch()

  return (
    <Box marginTop={-6}>
      <FormGroup row style={{alignItems: 'flex-end'}}>
        <Autocomplete
          options={searchOptions}
          getOptionLabel={(option) => option.name}
          style={{ width: 350 }}
          onChange={(_, value) => {
            if (value) {
              history.push(url + value.path)
            }
          }}
          renderInput={(params) => <TextField {...params} label="search" margin="normal" />}
        />
        <Box flexGrow={1} />
        <FormControlLabel
          control={
            <Checkbox
              checked={config.showCodeSpecific}
              onChange={handleConfigChange}
              name="showCodeSpecific"
            />
          }
          label="code specific"
        />
        <FormControlLabel
          control={
            <Checkbox
              checked={config.showAllDefined}
              onChange={handleConfigChange}
              name="showAllDefined"
            />
          }
          label="all defined"
        />
        <FormControlLabel
          control={
            <Checkbox
              checked={config.showMeta}
              onChange={handleConfigChange}
              name="showMeta" />
          }
          label="definitions"
        />
      </FormGroup>
    </Box>
  )
}
ArchiveConfigForm.propTypes = ({
  searchOptions: PropTypes.arrayOf(PropTypes.object).isRequired
})

function archiveAdaptorFactory(data, sectionDef) {
  return new SectionAdaptor(data, sectionDef || rootSections.find(def => def.name === 'EntryArchive'), {archive: data})
}

function archiveSearchOptions(data) {
  const options = []
  const optionDefs = {}
  const optionKeys = {}
  function traverse(data, def, parentPath) {
    for (let key in data) {
      const childDef = def._properties[key]
      if (!childDef) {
        continue
      }

      const child = data[key]
      if (!child) {
        continue
      }

      let path = `${parentPath}/${key}`
      if (childDef.m_def === 'SubSection') {
        const sectionDef = resolveRef(childDef.sub_section)
        if (Array.isArray(child) && child.length > 0 && child[0]) {
          if (child.length > 1) {
            child.forEach((value, index) => traverse(value, sectionDef, `${path}:${index}`))
          } else {
            traverse(child[0], sectionDef, path)
          }
        } else {
          traverse(child, sectionDef, path)
        }
      }

      if (optionDefs[childDef._qualifiedName]) {
        continue
      }
      optionDefs[childDef._qualifiedName] = childDef

      const option = {
        name: key,
        data: data,
        def: childDef,
        path: path
      }
      options.push(option)

      if (optionKeys[key]) {
        const addPath = option => {
          const parents = option.path.split('/')
          const parent = parents[parents.length - 2].replace(/:[0-9]+$/, '')
          option.name += ` (${parent})`
        }
        if (!optionKeys[key].name.includes('(')) {
          addPath(optionKeys[key])
        }
        addPath(option)
      } else {
        optionKeys[key] = option
      }
    }
  }
  traverse(data, rootSections.find(def => def.name === 'EntryArchive'), '')
  return options
}

class ArchiveAdaptor extends Adaptor {
  constructor(obj, def, context) {
    super(obj)
    this.def = def
    this.context = context
  }

  adaptorFactory(obj, def, context) {
    if (def.m_def === 'Section') {
      return new SectionAdaptor(obj, def, context || this.context)
    } else if (def.m_def === 'Quantity') {
      return new QuantityAdaptor(obj, def, context || this.context)
    }
  }

  itemAdaptor(key) {
    if (key === '_metainfo') {
      return metainfoAdaptorFactory(this.def)
    } else {
      throw new Error('Unknown item key')
    }
  }
}

class SectionAdaptor extends ArchiveAdaptor {
  itemAdaptor(key) {
    const [name, index] = key.split(':')
    const property = this.def._properties[name]
    const value = this.e[name]
    if (!property) {
      return super.itemAdaptor(key)
    } else if (property.m_def === 'SubSection') {
      const sectionDef = resolveRef(property.sub_section)
      if (property.repeats) {
        return this.adaptorFactory(value[parseInt(index || 0)], sectionDef)
      } else {
        return this.adaptorFactory(value, sectionDef)
      }
    } else if (property.m_def === 'Quantity') {
      if (property.type.type_kind === 'reference' && property.shape.length === 0) {
        return this.adaptorFactory(resolveRef(value, this.context.archive), resolveRef(property.type.type_data))
      }
      return this.adaptorFactory(value, property)
    } else {
      throw new Error('Unknown metainfo meta definition')
    }
  }
  render() {
    return <Section section={this.e} def={this.def} />
  }
}

class QuantityAdaptor extends ArchiveAdaptor {
  render() {
    return <Quantity value={this.e} def={this.def} />
  }
}

function QuantityItemPreview({value, def}) {
  if (def.type.type_kind === 'reference') {
    return <Box component="span" fontStyle="italic">
      <Typography component="span">reference ...</Typography>
    </Box>
  }
  if (def.shape.length > 0) {
    const dimensions = []
    let current = value
    for (let i = 0; i < def.shape.length; i++) {
      dimensions.push(current.length)
      current = current[0]
    }
    let typeLabel
    if (def.type.type_kind === 'python') {
      typeLabel = 'list'
    } else {
      if (dimensions.length === 1) {
        typeLabel = 'vector'
      } else if (dimensions.length === 2) {
        typeLabel = 'matrix'
      } else {
        typeLabel = 'tensor'
      }
    }
    return <Box component="span" whiteSpace="nowrap" fontStyle="italic">
      <Typography component="span">
        {dimensions.map((dimension, index) => (
          <span key={index}>
            {index > 0 && <span>&nbsp;&times;&nbsp;</span>}{String(dimension)}
          </span>
        ))}&nbsp;{typeLabel}
      </Typography>
    </Box>
  } else {
    return <Box component="span" whiteSpace="nowarp">
      <Number component="span" variant="body1" value={value} exp={8} />
      {def.unit && <Typography component="span">&nbsp;{def.unit}</Typography>}
    </Box>
  }
}
QuantityItemPreview.propTypes = ({
  value: PropTypes.any,
  def: PropTypes.object.isRequired
})

function QuantityValue({value, def}) {
  return <Box
    marginTop={2} marginBottom={2} textAlign="center" fontWeight="bold"
  >
    {def.shape.length > 0 ? <Matrix values={value} shape={def.shape} invert={def.shape.length === 1} /> : <Number value={value} exp={16} variant="body2" />}
    {def.shape.length > 0 &&
      <Typography noWrap variant="caption">
        ({def.shape.map((dimension, index) => <span key={index}>
          {index > 0 && <span>&nbsp;&times;&nbsp;</span>}{String(dimension)}
        </span>)}&nbsp;)
      </Typography>
    }
    {def.unit && <Typography noWrap>{def.unit}</Typography>}
  </Box>
}
QuantityValue.propTypes = ({
  value: PropTypes.any,
  def: PropTypes.object.isRequired
})

/**
 * Used to retrieve an optional visualization for a section displayed in the
 * browser.
 */
function getVisualization(section, def, archiveContext) {
  // Structure visualization for section_system
  if (def.name === 'section_system') {
    return <Structure viewer={archiveContext.structure} system={section}></Structure>
  }
}

function Section({section, def}) {
  const config = useRecoilValue(configState)
  const filter = config.showCodeSpecific ? def => true : def => !def.name.startsWith('x_')
  const archiveContext = useContext(ArchiveContext)

  return <Content>
    <Title def={def} data={section} kindLabel="section" />
    {getVisualization(section, def, archiveContext)}
    <Compartment title="sub sections">
      {def.sub_sections
        .filter(subSectionDef => section[subSectionDef.name] || config.showAllDefined)
        .filter(filter)
        .map(subSectionDef => {
          const key = subSectionDef.name
          const disabled = section[key] === undefined
          if (!disabled && subSectionDef.repeats && section[key].length > 1) {
            return <List
              key={subSectionDef.name}
              itemKey={subSectionDef.name}
              title={subSectionDef.name} disabled={disabled}
            />
          } else {
            return <Item key={key} itemKey={key} disabled={disabled}>
              <Typography component="span">
                <Box fontWeight="bold" component="span">
                  {subSectionDef.name}
                </Box>
              </Typography>
            </Item>
          }
        })
      }
    </Compartment>
    <Compartment title="quantities">
      {def.quantities
        .filter(quantityDef => section[quantityDef.name] || config.showAllDefined)
        .filter(filter)
        .map(quantityDef => {
          const key = quantityDef.name
          const disabled = section[key] === undefined
          return <Item key={key} itemKey={key} disabled={disabled}>
            <Box component="span" whiteSpace="nowrap">
              <Typography component="span">
                <Box fontWeight="bold" component="span">
                  {quantityDef.name}
                </Box>
              </Typography>{!disabled && <span>&nbsp;=&nbsp;<QuantityItemPreview value={section[quantityDef.name]} def={quantityDef} /></span>}
            </Box>
          </Item>
        })
      }
    </Compartment>
    <Meta def={def} />
  </Content>
}
Section.propTypes = ({
  section: PropTypes.object.isRequired,
  def: PropTypes.object.isRequired
})

function Quantity({value, def}) {
  return <Content>
    <Title def={def} data={value} kindLabel="value" />
    <QuantityValue value={value} def={def} />
    <Meta def={def} />
  </Content>
}
Quantity.propTypes = ({
  value: PropTypes.any,
  def: PropTypes.object.isRequired
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
  </Compartment>
}
Meta.propTypes = ({
  def: PropTypes.object
})
