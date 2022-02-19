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
import React, { useContext, useEffect, useMemo, useState } from 'react'
import PropTypes from 'prop-types'
import { atom, useRecoilState, useRecoilValue } from 'recoil'
import { Box, FormGroup, FormControlLabel, Checkbox, TextField, Typography, makeStyles, Tooltip, IconButton } from '@material-ui/core'
import { useRouteMatch, useHistory } from 'react-router-dom'
import Autocomplete from '@material-ui/lab/Autocomplete'
import Browser, { Item, Content, Compartment, Adaptor, formatSubSectionName, laneContext } from './Browser'
import { metainfoDef, QuantityMDef, resolveRef, rootSections, SectionMDef, SubSectionMDef } from './metainfo'
import { ArchiveTitle, metainfoAdaptorFactory, DefinitionLabel } from './MetainfoBrowser'
import { Matrix, Number } from './visualizations'
import Markdown from '../Markdown'
import { Overview } from './Overview'
import { toUnitSystem, useUnits } from '../../units'
import ArrowRightIcon from '@material-ui/icons/ArrowRight'
import ArrowDownIcon from '@material-ui/icons/ArrowDropDown'
import grey from '@material-ui/core/colors/grey'
import classNames from 'classnames'
import { useApi } from '../api'
import { useErrors } from '../errors'
import { SourceApiCall, SourceApiDialogButton } from '../buttons/SourceDialogButton'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import { Download } from '../entry/Download'
import SectionEditor from './SectionEditor'
import { useEntryContext } from '../entry/EntryContext'
import SaveIcon from '@material-ui/icons/Save'

export const configState = atom({
  key: 'config',
  default: {
    'showMeta': false,
    'showCodeSpecific': false,
    'showAllDefined': false
  }
})

const ArchiveBrowser = React.memo(({data}) => {
  const searchOptions = useMemo(() => archiveSearchOptions(data), [data])

  // For some reason, this hook does not work in all of the components used in
  // the Browser (notably: Quantity, QuantityItemPreview). In order to pass the
  // up-to-date unit information, we pass the hook value down the component
  // hierarchy.
  data.resources = data.resources || {}
  return (
    <Browser
      adaptor={archiveAdaptorFactory(data, undefined)}
      form={<ArchiveConfigForm searchOptions={searchOptions} data={data}/>}
    />
  )
})
ArchiveBrowser.propTypes = ({
  data: PropTypes.object.isRequired
})
export default ArchiveBrowser

export const ArchiveSaveButton = React.memo(function ArchiveSaveButton(props) {
  const {editable, archiveHasChanges, saveArchive} = useEntryContext()
  return <React.Fragment>
    {editable &&
      <IconButton
        disabled={!archiveHasChanges} color="primary"
        onClick={saveArchive}
      >
        <SaveIcon/>
      </IconButton>
    }
  </React.Fragment>
})

function ArchiveConfigForm({searchOptions, data}) {
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

  const entryId = data?.metadata?.entry_id

  return (
    <Box padding={0}>
      <FormGroup row style={{alignItems: 'center'}}>
        <Box style={{width: 350, height: 60}}>
          <Autocomplete
            options={searchOptions}
            getOptionLabel={(option) => option.name}
            style={{ width: 500, marginTop: -20 }}
            onChange={(_, value) => {
              if (value) {
                history.push(url + value.path)
              }
            }}
            renderInput={(params) => <TextField
              {...params} variant="filled"
              size="small" label="search" margin="normal"
            />}
          />
        </Box>
        <Box flexGrow={1} />
        <Tooltip title="Enable to also show all code specific data">
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
        </Tooltip>
        <Tooltip title="Enable to also show metadata that is in principle available, but not within this entry">
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
        </Tooltip>
        <Tooltip title="Show the Metainfo definition on the bottom of each lane">
          <FormControlLabel
            control={
              <Checkbox
                checked={config.showMeta}
                onChange={handleConfigChange}
                name="showMeta" />
            }
            label="definitions"
          />
        </Tooltip>
        {entryId && <Download
          tooltip="download the archive"
          url={`entries/${entryId}/archive/download`}
          fileName={`${entryId}.json`}
          component={IconButton}
        >
          <DownloadIcon />
        </Download>}
        <SourceApiDialogButton maxWidth="lg" fullWidth>
          <SourceApiCall />
        </SourceApiDialogButton>
        <ArchiveSaveButton/>
      </FormGroup>
    </Box>
  )
}
ArchiveConfigForm.propTypes = ({
  data: PropTypes.object.isRequired,
  searchOptions: PropTypes.arrayOf(PropTypes.object).isRequired
})

export function archiveAdaptorFactory(data, sectionDef) {
  return new SectionAdaptor(data, sectionDef || rootSections.find(def => def.name === 'EntryArchive'), undefined, {archive: data})
}

function archiveSearchOptions(data) {
  const options = []
  const optionDefs = {}
  function traverse(data, def, parentName, parentPath) {
    for (let key in data) {
      const childDef = def._properties[key]
      if (!childDef) {
        continue
      }

      const child = data[key]
      if (!child) {
        continue
      }
      const path = parentPath ? `${parentPath}/${key}` : `/${key}`
      const name = parentName ? `${parentName}.${childDef.name}` : childDef.name

      if (optionDefs[childDef._qualifiedName]) {
        continue
      }
      optionDefs[childDef._qualifiedName] = childDef
      const option = {
        name: name, // key
        data: data,
        def: childDef,
        path: path
      }
      options.push(option)

      if (childDef.m_def === SubSectionMDef) {
        const sectionDef = resolveRef(childDef.sub_section)
        if (Array.isArray(child) && child.length > 0 && child[0]) {
          if (child.length > 1) {
            child.forEach((value, index) => traverse(value, sectionDef, name, `${path}:${index}`))
          } else {
            traverse(child[0], sectionDef, name, path)
          }
        } else {
          traverse(child, sectionDef, name, path)
        }
      }
    }
  }
  traverse(data, rootSections.find(def => def.name === 'EntryArchive'), null, null)
  return options
}

class ArchiveAdaptor extends Adaptor {
  constructor(obj, def, parent, context) {
    super(obj)
    this.def = def
    this.parent = parent
    this.context = context
  }

  adaptorFactory(obj, def, parent, context) {
    if (def.m_def === SectionMDef) {
      if (obj.m_def) {
        // Override the def given by the schema with the potentially more specific
        // def given by the data
        def = metainfoDef(obj.m_def)
      }
      return new SectionAdaptor(obj, def, parent, context || this.context)
    } else if (def.m_def === QuantityMDef) {
      if (def.type.type_kind === 'reference') {
        return new ReferenceAdaptor(obj, def, parent, context || this.context)
      } else {
        return new QuantityAdaptor(obj, def, parent, context || this.context)
      }
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
    } else if (property.m_def === SubSectionMDef) {
      const sectionDef = resolveRef(property.sub_section)
      if (property.repeats) {
        return this.adaptorFactory(value[parseInt(index || 0)], sectionDef, this.e)
      } else {
        return this.adaptorFactory(value, sectionDef, this.e)
      }
    } else if (property.m_def === QuantityMDef) {
      // References: sections and quantities
      if (property.type.type_kind === 'reference') {
        let reference = null
        if (property.shape.length === 0) {
          reference = value
        } else if (property.shape.length === 1) {
          const indexStr = key.split(':')[1]
          const index = parseInt(indexStr)
          reference = value[index]
        }
        if (!reference) {
          return this.adaptorFactory(value, property, this.e)
        }
        const resolved = resolveRef(reference, this.context.archive)
        // some sections cannot be resolved, because they are not part of the archive
        // user_id->user is one example
        if (!resolved) {
          return this.adaptorFactory(reference, property, this.e)
        }
        const resolvedDef = resolveRef(property.type.type_data)
        return this.adaptorFactory(resolved, resolvedDef, this.e)
      }
      // Regular quantities
      return this.adaptorFactory(value, property, this.e)
    } else {
      throw new Error('Unknown metainfo meta definition')
    }
  }
  render() {
    return <Section section={this.e} def={this.def} parent={this.parent} />
  }
}

class ReferenceAdaptor extends ArchiveAdaptor {
  render() {
    return <Reference value={this.e} def={this.def} />
  }
}

class QuantityAdaptor extends ArchiveAdaptor {
  render() {
    return <Quantity value={this.e} def={this.def} />
  }
}

function QuantityItemPreview({value, def}) {
  const units = useUnits()
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
    const val = (def.type.type_data === 'nomad.metainfo.metainfo._Datetime' ? new Date(value).toLocaleString() : value)
    const [finalValue, finalUnit] = def.unit
      ? toUnitSystem(val, def.unit, units, true)
      : [val, def.unit]
    return <Box component="span" whiteSpace="nowarp">
      <Number component="span" variant="body1" value={finalValue} exp={8} />
      {finalUnit && <Typography component="span">&nbsp;{finalUnit}</Typography>}
    </Box>
  }
}
QuantityItemPreview.propTypes = ({
  value: PropTypes.any,
  def: PropTypes.object.isRequired
})

const QuantityValue = React.memo(function QuantityValue({value, def}) {
  const units = useUnits()
  const val = (def.type.type_data === 'nomad.metainfo.metainfo._Datetime' ? new Date(value).toLocaleString() : value)
  const [finalValue, finalUnit] = def.unit
    ? toUnitSystem(val, def.unit, units, true)
    : [val, def.unit]

  let isMathValue = def.type.type_kind === 'numpy'
  if (isMathValue) {
    if (def.shape.length > 0) {
      return <Box textAlign="center">
        <Matrix
          values={finalValue}
          shape={def.shape}
          invert={def.shape.length === 1}
          type={def.type.type_data}
        />
        <Typography noWrap variant="caption">
          ({def.shape.map((dimension, index) => <span key={index}>
            {index > 0 && <span>&nbsp;&times;&nbsp;</span>}{String(dimension)}
          </span>)}&nbsp;)
        </Typography>
        {finalUnit && <Typography noWrap>{finalUnit}</Typography>}
      </Box>
    } else {
      return <Number value={finalValue} exp={16} variant="body1" unit={finalUnit}/>
    }
  } else {
    if (Array.isArray(finalValue)) {
      return <Typography>
        <ul style={{margin: 0}}>
          {finalValue.map((value, index) =>
            <li key={index}>{typeof value === 'object' ? JSON.stringify(value) : value}</li>)}
        </ul>
      </Typography>
    } else {
      return <Typography>{finalValue}</Typography>
    }
  }
})
QuantityValue.propTypes = ({
  value: PropTypes.any,
  def: PropTypes.object.isRequired
})

function Section({section, def, parent}) {
  const config = useRecoilValue(configState)

  if (!section) {
    console.error('section is not available')
    return ''
  }

  const filter = config.showCodeSpecific ? def => true : def => !def.name.startsWith('x_')
  let sub_sections = def._allProperties.filter(prop => prop.m_def === SubSectionMDef)
  if (def.name === 'EntryArchive') {
    // put the most abstract data (last added data) first, e.g. results, metadata, workflow, run
    sub_sections = [...def.sub_sections]
    sub_sections.reverse()
  }
  const quantities = def._allProperties.filter(prop => prop.m_def === QuantityMDef)

  let contents
  if (def.name === 'Sample') {
    contents = <Compartment title="edit">
      <SectionEditor sectionDef={def} section={section} />
    </Compartment>
  } else {
    contents = <React.Fragment>
      <Compartment title="sub sections">
        {sub_sections
          .filter(subSectionDef => section[subSectionDef.name] || config.showAllDefined)
          .filter(filter)
          .map(subSectionDef => {
            const key = subSectionDef.name
            const disabled = section[key] === undefined
            if (!disabled && subSectionDef.repeats && section[key].length > 1) {
              return <SubSectionList
                key={subSectionDef.name}
                subSectionDef={subSectionDef}
                disabled={disabled}
              />
            } else {
              return <Item key={key} itemKey={key} disabled={disabled}>
                <Typography component="span">
                  <Box fontWeight="bold" component="span">
                    {formatSubSectionName(subSectionDef.name)}
                  </Box>
                </Typography>
              </Item>
            }
          })
        }
      </Compartment>
      <Compartment title="quantities">
        {quantities
          .filter(quantityDef => section[quantityDef.name] !== undefined || config.showAllDefined)
          .filter(filter)
          .map(quantityDef => {
            const key = quantityDef.name
            const disabled = section[key] === undefined
            if (!disabled && quantityDef.type.type_kind === 'reference' && quantityDef.shape.length === 1) {
              return <ReferenceValuesList key={key} quantityDef={quantityDef} />
            }
            return (
              <Item key={key} itemKey={key} disabled={disabled}>
                <Box component="span" whiteSpace="nowrap" style={{maxWidth: 100, overflow: 'ellipses'}}>
                  <Typography component="span">
                    <Box fontWeight="bold" component="span">
                      {quantityDef.name}
                    </Box>
                  </Typography>{!disabled &&
                    <span>&nbsp;=&nbsp;
                      <QuantityItemPreview
                        value={section[quantityDef.name]}
                        def={quantityDef}
                      />
                    </span>
                  }
                </Box>
              </Item>
            )
          })
        }
      </Compartment>
    </React.Fragment>
  }
  return <Content>
    <ArchiveTitle def={def} data={section} kindLabel="section" />
    <Overview section={section} def={def}/>
    {contents}
    <Meta def={def} />
  </Content>
}
Section.propTypes = ({
  section: PropTypes.object.isRequired,
  def: PropTypes.object.isRequired,
  parent: PropTypes.any
})

function SubSectionList({subSectionDef}) {
  const lane = useContext(laneContext)
  const label = useMemo(() => {
    let key = subSectionDef.more?.label_quantity
    if (!key) {
      const sectionDef = resolveRef(subSectionDef.sub_section)
      key = sectionDef.more?.label_quantity
      if (!key) {
        key = ['name', 'type', 'id'].find(key => (
          sectionDef._properties[key] && sectionDef._properties[key].m_def === QuantityMDef
        ))
      }
    }
    return item => {
      return key && item[key]
    }
  }, [subSectionDef])
  const values = useMemo(() => lane.adaptor.e[subSectionDef.name].map(label), [lane.adaptor.e, subSectionDef.name, label])
  return <PropertyValuesList
    values={values}
    label={formatSubSectionName(subSectionDef.name) || 'list'} />
}
SubSectionList.propTypes = ({
  subSectionDef: PropTypes.object.isRequired
})

function ReferenceValuesList({quantityDef}) {
  const lane = useContext(laneContext)
  const values = useMemo(() => lane.adaptor.e[quantityDef.name].map(() => null), [lane.adaptor.e, quantityDef.name])
  return <PropertyValuesList
    values={values}
    label={quantityDef.name} />
}
ReferenceValuesList.propTypes = ({
  quantityDef: PropTypes.object.isRequired
})

const usePropertyValuesListStyles = makeStyles(theme => ({
  title: {
    color: theme.palette.text.primary,
    textDecoration: 'none',
    margin: `0 -${theme.spacing(1)}px`,
    whiteSpace: 'nowrap',
    display: 'flex',
    fontWeight: 'bold'
  },
  selected: {
    backgroundColor: theme.palette.primary.main,
    color: theme.palette.primary.contrastText,
    whiteSpace: 'nowrap'
  },
  unSelected: {
    '&:hover': {
      backgroundColor: grey[300]
    }
  }
}))
function PropertyValuesList({label, values}) {
  const classes = usePropertyValuesListStyles()
  const [open, setOpen] = useState(false)
  const lane = useContext(laneContext)
  const selected = lane.next && lane.next.key

  return <div>
    <Typography onClick={() => setOpen(!open)} className={classNames(
      classes.title,
      (!open && selected && selected.startsWith(label + ':')) ? classes.selected : classes.unSelected
    )}>
      {open ? <ArrowDownIcon/> : <ArrowRightIcon/>}
      <span>{label}</span>
    </Typography>
    {open &&
      <div>
        {values.map((item, index) => (
          <Item key={index} itemKey={`${label}:${index}`}>
            <Box component="span" marginLeft={2}>
              <Typography component="span">{item || index}</Typography>
            </Box>
          </Item>
        ))}
      </div>
    }
  </div>
}
PropertyValuesList.propTypes = ({
  label: PropTypes.string.isRequired,
  values: PropTypes.arrayOf(PropTypes.string).isRequired
})

function Quantity({value, def}) {
  return <Content>
    <ArchiveTitle def={def} data={value} kindLabel="value" />
    <Compartment title="value">
      <QuantityValue
        value={value}
        def={def}
      />
    </Compartment>
    <Meta def={def} />
  </Content>
}
Quantity.propTypes = ({
  value: PropTypes.any,
  def: PropTypes.object.isRequired
})

function Reference({value, def}) {
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [loading, setLoading] = useState(true)
  const {data, update} = useContext(laneContext)
  useEffect(() => {
    const url = value.split('#')[0]
    const upload_id = data.metadata.upload_id
    if (data.resources[url]) {
      setLoading(false)
      return
    }

    if (!(url.startsWith('../upload/archive/') && upload_id)) {
      setLoading(false)
      return
    }

    api.get(`uploads/${upload_id}/${url.slice('../upload/'.length)}`)
      .then(response => {
        data.resources[url] = response.data.archive
        update()
      })
      .catch(raiseError)
  }, [api, data.metadata.upload_id, data.resources, raiseError, setLoading, update, value])

  if (loading) {
    return <Content>
      <Typography>loading ...</Typography>
    </Content>
  }

  return <Content>
    <ArchiveTitle def={def} data={value} kindLabel="value" />
    <Compartment title="reference">
      <Typography color="error">Cannot resolve reference.</Typography>
      <Typography>{value}</Typography>
    </Compartment>
    <Meta def={def} />
  </Content>
}
Reference.propTypes = ({
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
