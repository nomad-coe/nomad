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
import React, { useCallback, useContext, useEffect, useMemo, useState } from 'react'
import PropTypes from 'prop-types'
import { atom, useRecoilState, useRecoilValue } from 'recoil'
import { Box, FormGroup, FormControlLabel, Checkbox, TextField, Typography, makeStyles, Tooltip, IconButton, useTheme, Grid } from '@material-ui/core'
import { useRouteMatch, useHistory } from 'react-router-dom'
import Autocomplete from '@material-ui/lab/Autocomplete'
import Browser, { Item, Content, Compartment, Adaptor, formatSubSectionName, laneContext, useLane } from './Browser'
import { RawFileAdaptor } from './FileBrowser'
import { isEditable, metainfoDef, QuantityMDef, removeSubSection, resolveRef, rootSections, SectionMDef, SubSectionMDef } from './metainfo'
import { ArchiveTitle, metainfoAdaptorFactory, DefinitionLabel } from './MetainfoBrowser'
import { Matrix, Number } from './visualizations'
import Markdown from '../Markdown'
import { Overview } from './Overview'
import { Quantity as Q, useUnits } from '../../units'
import ArrowRightIcon from '@material-ui/icons/ArrowRight'
import ArrowDownIcon from '@material-ui/icons/ArrowDropDown'
import grey from '@material-ui/core/colors/grey'
import classNames from 'classnames'
import { useApi } from '../api'
import { useErrors } from '../errors'
import { SourceApiCall, SourceApiDialogButton, SourceJsonDialogButton } from '../buttons/SourceDialogButton'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import { Download } from '../entry/Download'
import SectionEditor from './SectionEditor'
import { useEntryContext } from '../entry/EntryContext'
import SaveIcon from '@material-ui/icons/Save'
import AddIcon from '@material-ui/icons/AddCircle'
import CodeIcon from '@material-ui/icons/Code'
import DeleteIcon from '@material-ui/icons/Delete'
import { getLineStyles } from '../../utils'
import Plot from '../visualization/Plot'

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
          url={`entries/${entryId}/archive/download?ignore_mime_type=true`}
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
    context = context || this.context
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
  async itemAdaptor(key, api) {
    const [name, index] = key.split(':')
    const property = this.def._properties[name]
    const value = this.e[name]
    if (!property) {
      return super.itemAdaptor(key)
    } else if (property.m_def === SubSectionMDef) {
      const sectionDef = resolveRef(property.sub_section)
      let subSectionAdaptor
      let subSectionIndex = -1
      if (property.repeats) {
        subSectionIndex = parseInt(index || 0)
        subSectionAdaptor = this.adaptorFactory(value[subSectionIndex], sectionDef, this.e)
      } else {
        subSectionAdaptor = this.adaptorFactory(value, sectionDef, this.e)
      }
      subSectionAdaptor.parentRelation = {
        parent: this.e,
        subSectionDef: property,
        subSectionIndex: subSectionIndex
      }
      return subSectionAdaptor
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
        const context = {
          ...this.context,
          isReferenced: true
        }
        if (reference.includes('#')) {
          const [url] = reference.split('#')
          if (url !== '') {
            context.archive = this.context.archive.resources[url]
          }
        }
        return this.adaptorFactory(resolved, resolvedDef, this.e, context)
      }
      // Regular quantities
      if (property.m_annotations?.browser) {
        if (property.m_annotations.browser[0].adaptor === 'RawFileAdaptor') {
          const uploadId = this.context.archive.metadata.upload_id
          const path = this.e[property.name]
          const response = await api.get(`uploads/${uploadId}/rawdir/${path}`)
          return new RawFileAdaptor(uploadId, path, response.file_metadata, false)
        }
      }
      return this.adaptorFactory(value, property, this.e)
    } else {
      throw new Error('Unknown metainfo meta definition')
    }
  }
  render() {
    return <Section
      section={this.e}
      def={this.def}
      parent={this.parent} parentRelation={this.parentRelation} />
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
  if (def.m_annotations?.eln?.[0]?.component === 'RichTextEditQuantity') {
    return <Box component="span" whiteSpace="nowrap" fontStyle="italic">
      <Typography component="span">rich text</Typography>
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
    let finalValue = (def.type.type_data === 'nomad.metainfo.metainfo._Datetime' ? new Date(value).toLocaleString() : value)
    let finalUnit
    if (def.unit) {
      const a = new Q(finalValue, def.unit).toSystem(units)
      finalValue = a.value
      finalUnit = a.unit.label
    }
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

  const getRenderValue = useCallback(value => {
    let finalValue = (def.type.type_data === 'nomad.metainfo.metainfo._Datetime' ? new Date(value).toLocaleString() : value)
    let finalUnit
    if (def.unit) {
      const a = new Q(finalValue, def.unit).toSystem(units)
      finalValue = a.value
      finalUnit = a.unit.label
    }
    return [finalValue, finalUnit]
  }, [def, units])

  const isMathValue = def.type.type_kind === 'numpy'
  if (isMathValue) {
    const [finalValue, finalUnit] = getRenderValue(value)
    if (def.shape.length > 0) {
      return <Box textAlign="center">
        <Matrix
          values={finalValue}
          shape={def.shape}
          invert={def.shape.length === 1}
          type={def.type.type_data}
          key={`matrix:${def.name}`}
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
  } else if (def.m_annotations?.eln?.[0]?.component === 'RichTextEditQuantity') {
    return <div dangerouslySetInnerHTML={{ __html: value }}/>
  } else {
    if (Array.isArray(value)) {
      return <ul style={{margin: 0}}>
        {value.map((value, index) => {
          const [finalValue] = getRenderValue(value)
          return <li key={index}>
            <Typography>{typeof finalValue === 'object' ? JSON.stringify(finalValue) : finalValue?.toString()}</Typography>
          </li>
        })}
      </ul>
    } else {
      const [finalValue] = getRenderValue(value)
      return <Typography>{typeof finalValue === 'object' ? JSON.stringify(finalValue) : finalValue?.toString()}</Typography>
    }
  }
})
QuantityValue.propTypes = ({
  value: PropTypes.any,
  def: PropTypes.object.isRequired
})

function Section({section, def, parentRelation}) {
  const {editable, handleArchiveChanged} = useEntryContext() || {}
  const config = useRecoilValue(configState)
  const [showJson, setShowJson] = useState(false)
  const lane = useLane()
  const history = useHistory()

  const sectionIsEditable = useMemo(() => {
    return editable && isEditable(def) && !lane.adaptor.context.isReferenced
  }, [editable, def, lane])

  const actions = useMemo(() => {
    if (!sectionIsEditable) {
      return <SourceJsonDialogButton
        buttonProps={{size: 'small'}}
        tooltip={`Show section data as JSON`}
        title={`Underlying section data as JSON`}
        data={section}
      />
    }

    const handleDelete = () => {
      removeSubSection(
        parentRelation.parent,
        parentRelation.subSectionDef,
        parentRelation.subSectionIndex)
      handleArchiveChanged()
      history.push(lane.prev.path)
    }

    return <Grid container justifyContent="space-between" wrap="nowrap" spacing={1}>
      <Grid item>
        <IconButton onClick={() => setShowJson(value => !value)} size="small">
          <CodeIcon />
        </IconButton>
      </Grid>
      <Grid item>
        <IconButton onClick={handleDelete} size="small">
          <DeleteIcon />
        </IconButton>
      </Grid>
    </Grid>
  }, [setShowJson, sectionIsEditable, parentRelation, lane, history, handleArchiveChanged, section])

  const renderQuantity = useCallback(quantityDef => {
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
  }, [section])

  if (!section) {
    console.error('section is not available')
    return ''
  }

  const filter = config.showCodeSpecific ? def => !def.virtual : def => !def.virtual && !def.name.startsWith('x_')
  let sub_sections = def._allProperties.filter(prop => prop.m_def === SubSectionMDef)
  if (def.name === 'EntryArchive') {
    // put the most abstract data (last added data) first, e.g. results, metadata, workflow, run
    sub_sections = [...def.sub_sections]
    sub_sections.reverse()
  }
  const quantities = def._allProperties.filter(prop => prop.m_def === QuantityMDef)

  const subSectionCompartment = (
    <Compartment title="sub sections">
      {sub_sections
        .filter(subSectionDef => section[subSectionDef.name] || config.showAllDefined || sectionIsEditable)
        .filter(filter)
        .map(subSectionDef => {
          return <SubSection
            key={subSectionDef.name}
            subSectionDef={subSectionDef}
            section={section}
            editable={sectionIsEditable}
          />
        })
      }
    </Compartment>
  )

  let contents
  if (sectionIsEditable) {
    contents = <React.Fragment>
      {quantities.length > 0 && (
        <Compartment title="quantities">
          <SectionEditor sectionDef={def} section={section} showJson={showJson} />
          <Box marginTop={2}>
            {quantities
              .filter(filter)
              .filter(quantityDef => !quantityDef.m_annotations?.eln)
              .map(renderQuantity)
            }
          </Box>
        </Compartment>
      )}
      {subSectionCompartment}
    </React.Fragment>
  } else {
    contents = <React.Fragment>
      {subSectionCompartment}
      <Compartment title="quantities">
        {quantities
          .filter(quantityDef => section[quantityDef.name] !== undefined || config.showAllDefined)
          .filter(filter)
          .map(renderQuantity)
        }
      </Compartment>
    </React.Fragment>
  }
  const eln = def?.m_annotations?.eln
  const laneWidth = (eln && eln.length > 0 ? eln[0].lane_width : undefined)
  const otherProps = (laneWidth ? {minWidth: laneWidth, maxWidth: laneWidth} : undefined)
  return <Content {...otherProps}>
    <ArchiveTitle def={def} data={section} kindLabel="section" actions={actions} />
    <Overview section={section} def={def}/>
    {contents}
    <Meta def={def} />
  </Content>
}
Section.propTypes = ({
  section: PropTypes.object.isRequired,
  def: PropTypes.object.isRequired,
  subSection: PropTypes.object,
  parentRelation: PropTypes.object
})

function SubSection({subSectionDef, section, editable}) {
  const {handleArchiveChanged} = useEntryContext() || {}
  const lane = useLane()
  const history = useHistory()
  const { label, getItemLabel } = useMemo(() => {
    const sectionDef = resolveRef(subSectionDef.sub_section)
    let itemLabelKey = sectionDef.more?.label_quantity
    if (!itemLabelKey) {
      itemLabelKey = ['name', 'type', 'id'].find(key => (
        sectionDef._properties[key] && sectionDef._properties[key].m_def === QuantityMDef
      ))
    }
    const labelQuantity = itemLabelKey && sectionDef._properties[itemLabelKey]
    const getItemLabel = item => {
      if (labelQuantity) {
        const value = item[itemLabelKey]
        if (value) {
          return <QuantityValue value={item[itemLabelKey]} def={labelQuantity}/>
        }
      } else if (itemLabelKey) {
        return item[itemLabelKey]
      }
      return null
    }
    return {
      label: formatSubSectionName(subSectionDef.name),
      getItemLabel: getItemLabel
    }
  }, [subSectionDef])

  const handleAdd = useCallback(() => {
    let subSectionKey = subSectionDef.name
    if (subSectionDef.repeats) {
      let values = section[subSectionDef.name]
      if (!values) {
        values = []
        section[subSectionDef.name] = values
      }
      values.push({})
      if (values.length > 1) {
        subSectionKey += `:${values.length - 1}`
      }
    } else {
      section[subSectionDef.name] = {}
    }
    handleArchiveChanged()
    history.push(`${lane.path}/${subSectionKey}`)
  }, [subSectionDef, section, lane, history, handleArchiveChanged])

  const values = section[subSectionDef.name]
  const showList = subSectionDef.repeats && values && values.length > 1
  const actions = editable && (subSectionDef.repeats || !values) && (
    <Box marginRight={!showList && values ? -1 : 2}>
      <IconButton onClick={handleAdd} size="small">
        <AddIcon style={{fontSize: 20}} />
      </IconButton>
    </Box>
  )

  if (showList) {
    return <PropertyValuesList
      label={label || 'list'} actions={actions}
      values={(section[subSectionDef.name] || []).map(getItemLabel)}
    />
  } else {
    return (
      <Item
        itemKey={subSectionDef.name} disabled={!values}
        actions={actions}
      >
        <Typography component="span">
          <Box fontWeight="bold" component="span">
            {label}
          </Box>
        </Typography>
      </Item>
    )
  }
}
SubSection.propTypes = ({
  subSectionDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  editable: PropTypes.bool
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
  root: {
    margin: `0 -${theme.spacing(1)}px`,
    padding: `0 ${theme.spacing(1)}px`,
    width: '100%',
    display: 'flex',
    flexDirection: 'row',
    wrap: 'nowrap',
    alignItems: 'center'
  },
  title: {
    flexGrow: 1,
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
  },
  actions: {}
}))
function PropertyValuesList({label, values, actions}) {
  const classes = usePropertyValuesListStyles()
  const [open, setOpen] = useState(false)
  const lane = useContext(laneContext)
  const selected = lane.next && lane.next.key
  const showSelected = !open && selected && selected.startsWith(label + ':')
  return <React.Fragment>
    <div className={classNames(
      classes.root, showSelected ? classes.selected : classes.unSelected)}
    >
      <Typography onClick={() => setOpen(!open)} className={classes.title}>
        {open ? <ArrowDownIcon/> : <ArrowRightIcon/>}
        <span role="item-list">{label}</span>
      </Typography>
      {actions && <div className={classes.actions}>
        {actions}
      </div>}
    </div>
    {open &&
      <div data-testid={`item-list:${label}`} >
        {values.map((item, index) => (
          <Item key={index} itemKey={`${label}:${index}`}>
            <Box display="flex" flexDirection="row" flexGrow={1}>
              <Box component="span" marginLeft={2}>
                { item && typeof item === 'object'
                  ? item // item should be a react component
                  : <Typography component="span">{item || index}</Typography>
                }
              </Box>
            </Box>
          </Item>
        ))}
      </div>
    }
  </React.Fragment>
}
PropertyValuesList.propTypes = ({
  label: PropTypes.string.isRequired,
  values: PropTypes.arrayOf(PropTypes.object).isRequired,
  onAdd: PropTypes.func,
  onRemove: PropTypes.func,
  actions: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
})

const XYPlot = React.memo(function XYPlot({plot, section, sectionDef, title}) {
  const theme = useTheme()
  const units = useUnits()

  const [data, layout] = useMemo(() => {
    const toUnit = quantityDef => {
      const value = section[quantityDef.name]
      const unit = quantityDef.unit
      if (unit) {
        const quantity = new Q(value, quantityDef.unit).toSystem(units)
        return [quantity.value, quantity.unit.label]
      } else {
        return [value, unit]
      }
    }
    const [xValues, xUnit] = toUnit(sectionDef._properties[plot.x])
    const [yValues, yUnit] = toUnit(sectionDef._properties[plot.y])

    const data = [
      {
        x: xValues,
        y: yValues,
        type: 'scatter',
        mode: 'lines',
        line: getLineStyles(1, theme)
      }
    ]

    const layout = {
      yaxis: {
        title: {
          text: yUnit ? `${plot.y} (${yUnit})` : plot.y
        }
      },
      xaxis: {
        title: {
          text: yUnit ? `${plot.x} (${xUnit})` : plot.x
        }
      }
    }

    return [data, layout]
  }, [section, plot, sectionDef, theme, units])

  return <Box minWidth={500}>
    <Plot
      data={data}
      layout={layout}
      floatTitle={title}
      fixedMargins={true}
    />
  </Box>
})
XYPlot.propTypes = {
  plot: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  sectionDef: PropTypes.object.isRequired,
  title: PropTypes.string
}

function Quantity({value, def}) {
  const {prev} = useLane()
  return <Content>
    <ArchiveTitle def={def} data={value} kindLabel="value" />
    {def.m_annotations?.plot && (
      <Compartment title="plot">
        <XYPlot
          section={prev.adaptor.e}
          sectionDef={prev.adaptor.def}
          plot={def.m_annotations?.plot[0]}
          title={def.name}
        />
      </Compartment>
    )}
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
  const {adaptor, update} = useContext(laneContext)
  const archive = adaptor.context.archive
  useEffect(() => {
    const url = value.split('#')[0]
    const upload_id = archive.metadata.upload_id
    if (archive.resources[url]) {
      setLoading(false)
      return
    }

    if (!(url.startsWith('../upload/archive/') && upload_id)) {
      setLoading(false)
      return
    }

    api.get(`uploads/${upload_id}/${url.slice('../upload/'.length)}`)
      .then(response => {
        archive.resources[url] = response.data.archive
        update()
      })
      .catch(raiseError)
  }, [api, archive.metadata.upload_id, archive.resources, raiseError, setLoading, update, value])

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
