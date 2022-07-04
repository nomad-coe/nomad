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
import {
  Box, FormGroup, FormControlLabel, Checkbox, TextField, Typography, makeStyles, Tooltip,
  IconButton, useTheme, Grid, Dialog, DialogContent, DialogContentText, DialogActions,
  Button, MenuItem } from '@material-ui/core'
import {useRouteMatch, useHistory} from 'react-router-dom'
import Autocomplete from '@material-ui/lab/Autocomplete'
import Browser, { Item, Content, Compartment, Adaptor, formatSubSectionName, laneContext, useLane } from './Browser'
import { RawFileAdaptor } from './FileBrowser'
import {
  isEditable, PackageMDef, QuantityMDef, removeSubSection, resolveRef, resolveRefAsync, SectionMDef, SubSectionMDef,
  useMetainfo
} from './metainfo'
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
import {getLineStyles, titleCase} from '../../utils'
import Plot from '../visualization/Plot'
import { useUploadPageContext } from '../uploads/UploadPageContext'
import {EntryButton} from '../nav/Routes'
import NavigateIcon from '@material-ui/icons/MoreHoriz'
import {ErrorHandler} from '../ErrorHandler'
import Alert from '@material-ui/lab/Alert'
import _ from 'lodash'
import ReloadIcon from '@material-ui/icons/Replay'
import UploadIcon from '@material-ui/icons/CloudUpload'

export function useBrowserAdaptorContext(data) {
  const entryContext = useEntryContext()
  const uploadPageContext = useUploadPageContext()
  const metainfo = useMetainfo(data)
  const {api} = useApi()

  const entryId = entryContext?.entryId
  const uploadId = entryContext?.uploadId || uploadPageContext?.uploadId
  const mainfile = entryContext?.metadata?.mainfile

  const context = useMemo(() => ({
    api: api,
    metainfo: metainfo,
    archive: data,
    resources: {},
    uploadId: uploadId,
    entryId: entryId,
    mainfile: mainfile
  }), [uploadId, entryId, mainfile, api, data, metainfo])

  return context
}

export const configState = atom({
  key: 'config',
  default: {
    'showMeta': false,
    'showCodeSpecific': false,
    'showAllDefined': false
  }
})

const ArchiveBrowser = React.memo(({data}) => {
  const context = useBrowserAdaptorContext(data)
  const {metainfo} = context
  const searchOptions = useMemo(() => {
    return metainfo ? archiveSearchOptions(data, metainfo) : []
  }, [data, metainfo])

  const adaptor = useMemo(() => {
    if (!context.metainfo) {
      return null
    }
    return archiveAdaptorFactory(context, data, undefined)
  }, [context, data])

  if (!adaptor) {
    return ''
  }

  // For some reason, this hook does not work in all of the components used in
  // the Browser (notably: Quantity, QuantityItemPreview). In order to pass the
  // up-to-date unit information, we pass the hook value down the component
  // hierarchy.
  context.resources = context.resources || {}
  context.archive = data
  return (
    <Browser
      adaptor={adaptor}
      form={<ArchiveConfigForm searchOptions={searchOptions} data={data}/>}
    />
  )
})
ArchiveBrowser.propTypes = ({
  data: PropTypes.object.isRequired
})
export default ArchiveBrowser

export const ArchiveSaveButton = React.memo(function ArchiveSaveButton(props) {
  const {editable, archiveHasChanges, saveArchive, reload} = useEntryContext()
  const [openErrorDialog, setOpenErrorDialog] = useState(false)
  const [disabled, setDisabled] = useState(false)

  const handleClick = useCallback(() => {
    saveArchive().catch(error => {
      if (error?.status === 409) {
        setOpenErrorDialog(true)
        setDisabled(true)
      }
    })
  }, [saveArchive])

  const handleReload = useCallback(() => {
    reload()
    setOpenErrorDialog(false)
  }, [reload])

  return <React.Fragment>
    {editable &&
      <IconButton
        disabled={!archiveHasChanges || disabled} color="primary"
        onClick={handleClick}
      >
        <Tooltip title="Save archive">
          <SaveIcon/>
        </Tooltip>
      </IconButton>
    }
    <Dialog
      open={openErrorDialog}
      aria-describedby="alert-dialog-description"
    >
      <DialogContent>
        <DialogContentText>
          The changes cannot be saved. The content has been modified by someone else.
        </DialogContentText>
      </DialogContent>
      <DialogActions>
        <Button onClick={() => setOpenErrorDialog(false)}>OK</Button>
        <Button onClick={handleReload} autoFocus>Reload</Button>
      </DialogActions>
    </Dialog>
  </React.Fragment>
})

export const ArchiveReloadButton = React.memo(function ArchiveSaveButton(props) {
  const {reload} = useEntryContext()

  return <React.Fragment>
    <IconButton
      color="primary"
      onClick={reload}
    >
      <Tooltip title="Reload archive">
        <ReloadIcon/>
      </Tooltip>
    </IconButton>
  </React.Fragment>
})

export const ArchiveDeleteButton = React.memo(function ArchiveDeleteButton(props) {
  const history = useHistory()
  const {editable, uploadId, entryId} = useEntryContext()
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [openDeleteConfirmDialog, setOpenDeleteConfirmDialog] = useState(false)

  const handleClick = useCallback(() => {
    setOpenDeleteConfirmDialog(true)
  }, [setOpenDeleteConfirmDialog])

  const handleDelete = useCallback(includeParentFolders => {
    setOpenDeleteConfirmDialog(false)
    const requestBody = {query: {entry_id: entryId}, include_parent_folders: includeParentFolders}
    api.post(`uploads/${uploadId}/action/delete-entry-files`, requestBody)
      .then(results => {
        history.push(`/user/uploads/upload/id/${uploadId}`)
      })
      .catch(err =>
        raiseError(err)
      )
  }, [uploadId, entryId, history, api, raiseError, setOpenDeleteConfirmDialog])

  return editable ? (
    <React.Fragment>
      <IconButton color="primary" onClick={handleClick}>
        <Tooltip title="Delete archive">
          <DeleteIcon/>
        </Tooltip>
      </IconButton>
      <Dialog
        open={openDeleteConfirmDialog}
        aria-describedby="alert-dialog-description"
      >
        <DialogContent>
          <DialogContentText id="alert-dialog-description">
            <b>Really delete this entry?</b>
          </DialogContentText>
          <DialogContentText>
            You can choose to delete only the mainfile, or to delete the mainfile and its folder.
          </DialogContentText>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setOpenDeleteConfirmDialog(false)} autoFocus>Cancel</Button>
          <Button onClick={() => handleDelete(false)}>Delete mainfile</Button>
          <Button onClick={() => handleDelete(true)}>Delete mainfile and folder</Button>
        </DialogActions>
      </Dialog>
    </React.Fragment>) : ''
})

const ArchiveConfigForm = React.memo(function ArchiveConfigForm({searchOptions, data}) {
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
        <ArchiveReloadButton />
        <ArchiveSaveButton/>
      </FormGroup>
    </Box>
  )
})
ArchiveConfigForm.propTypes = ({
  data: PropTypes.object.isRequired,
  searchOptions: PropTypes.arrayOf(PropTypes.object).isRequired
})

export const ArchiveReUploadButton = React.memo((props) => {
  const {uploadId, metadata, reload} = useEntryContext()
  const {api} = useApi()
  const {raiseError} = useErrors()

  const handleClick = useCallback((files) => {
    const input = document.createElement('input')
    input.type = 'file'
    input.onchange = (event) => {
      const file = event.target.files[0]
      if (!file) {
        return
      }
      const formData = new FormData() // eslint-disable-line no-undef
      formData.set('file', file, metadata.entry_name)
      api.put(`uploads/${uploadId}/raw/?wait_for_processing=true`, formData)
        .then(() => {
          reload()
        })
        .catch(raiseError)
    }

    input.click()
  }, [api, metadata.entry_name, raiseError, reload, uploadId])

  return <IconButton onClick={handleClick}>
    <Tooltip title="Replace this entry's mainfile">
      <UploadIcon/>
    </Tooltip>
  </IconButton>
})

export function archiveAdaptorFactory(context, data, sectionDef) {
  return new SectionAdaptor(
    {archive: data, ...context},
    data,
    sectionDef || context.metainfo?.getEntryArchiveDefinition(),
    undefined)
}

function archiveSearchOptions(data, metainfo) {
  const options = []
  const optionDefs = {}
  function traverse(data, def, parentName, parentPath) {
    for (const key in data) {
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
        const sectionDef = childDef.sub_section
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
  traverse(data, metainfo.getEntryArchiveDefinition(), null, null)
  return options
}

class ArchiveAdaptor extends Adaptor {
  constructor(context, obj, def, parent) {
    super(context)
    this.obj = obj
    this.def = def
    this.parent = parent

    if (!this.def) {
      throw new Error('Definitions must be given.')
    }
  }

  async adaptorFactory(obj, def, parent, childContext) {
    const context = childContext || this.context
    if (def === await this.context.metainfo.resolveDefinition(PackageMDef)) {
      return metainfoAdaptorFactory(this.context, obj)
    }

    if (def.m_def === SectionMDef) {
      if (obj.m_def) {
        // Override the def given by the schema with the potentially more specific
        // def given by the data
        def = await context.metainfo.resolveDefinition(obj.m_def, context)
      }
      return new SectionAdaptor(context, obj, def, parent)
    }

    if (def.m_def === QuantityMDef) {
      if (def.type.type_kind === 'reference') {
        return new ReferenceAdaptor(context, obj, def, parent)
      }

      return new QuantityAdaptor(context, obj, def, parent)
    }

    throw new Error('not implemented')
  }

  async itemAdaptor(key) {
    if (key === '_metainfo') {
      return metainfoAdaptorFactory(this.context, this.def)
    } else {
      throw new Error('Unknown item key')
    }
  }
}

class SectionAdaptor extends ArchiveAdaptor {
  async itemAdaptor(key) {
    const [name, index] = key.split(':')
    const property = this.def._properties[name]
    const value = this.obj[name]
    if (!property) {
      return super.itemAdaptor(key)
    } else if (property.m_def === SubSectionMDef) {
      const sectionDef = property.sub_section
      let subSectionAdaptor
      let subSectionIndex = -1
      if (property.repeats) {
        subSectionIndex = parseInt(index || 0)
        subSectionAdaptor = await this.adaptorFactory(value[subSectionIndex], sectionDef, this.obj)
      } else {
        subSectionAdaptor = await this.adaptorFactory(value, sectionDef, this.obj)
      }
      subSectionAdaptor.parentRelation = {
        parent: this.obj,
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
          return this.adaptorFactory(value, property, this.obj)
        }
        const childContext = {
          ...this.context,
          isReferenced: true
        }
        const resolved = await resolveRefAsync(reference, this.context.archive, this.context, archive => {
          childContext.archive = archive
        })
        // some sections cannot be resolved, because they are not part of the archive
        // user_id->user is one example
        if (!resolved) {
          return this.adaptorFactory(reference, property, this.obj)
        }
        const resolvedDef = property.type._referencedSection
        return this.adaptorFactory(resolved, resolvedDef, this.obj, childContext)
      }
      // Regular quantities
      if (property.m_annotations?.browser) {
        if (property.m_annotations.browser[0].adaptor === 'RawFileAdaptor') {
          const uploadId = this.context.archive.metadata.upload_id
          const path = this.obj[property.name]
          const response = await this.context.api.get(`uploads/${uploadId}/rawdir/${path}`)
          return new RawFileAdaptor(this.context, uploadId, path, response.file_metadata, false)
        }
      }
      return this.adaptorFactory(value, property, this.obj)
    } else {
      throw new Error('Unknown metainfo meta definition')
    }
  }
  render() {
    return <Section
      section={this.obj}
      def={this.def}
      parent={this.parent} parentRelation={this.parentRelation} />
  }
}

class ReferenceAdaptor extends ArchiveAdaptor {
  render() {
    return <Reference value={this.obj} def={this.def} />
  }
}

class QuantityAdaptor extends ArchiveAdaptor {
  render() {
    return <Quantity value={this.obj} def={this.def} />
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
      finalValue = a.value()
      finalUnit = a.label()
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
      finalValue = a.value()
      finalUnit = a.label()
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

  const navEntryId = useMemo(() => {
    return lane?.adaptor?.context?.archive?.metadata?.entry_id
  }, [lane])

  const sectionIsEditable = useMemo(() => {
    return editable && isEditable(def) && !lane.adaptor.context.isReferenced
  }, [editable, def, lane])

  const actions = useMemo(() => {
    const navButton = navEntryId && (
      <Grid item>
        <EntryButton entryId={navEntryId} component={IconButton} size="small">
          <NavigateIcon />
        </EntryButton>
      </Grid>
    )

    const jsonButton = !sectionIsEditable ? (
      <Grid item>
        <SourceJsonDialogButton
          buttonProps={{size: 'small'}}
          tooltip={`Show section data as JSON`}
          title={`Underlying section data as JSON`}
          data={section}
        />
      </Grid>
    ) : (
      <Grid item>
        <IconButton onClick={() => setShowJson(value => !value)} size="small">
          <CodeIcon />
        </IconButton>
      </Grid>
    )

    const handleDelete = () => {
      removeSubSection(
        parentRelation.parent,
        parentRelation.subSectionDef,
        parentRelation.subSectionIndex)
      handleArchiveChanged()
      history.push(lane.prev.path)
    }

    const deleteButton = sectionIsEditable && (
      <Grid item>
        <IconButton onClick={handleDelete} size="small">
          <DeleteIcon />
        </IconButton>
      </Grid>
    )

    return <Grid container justifyContent="space-between" wrap="nowrap" spacing={1}>
      {navButton}{jsonButton}{deleteButton}
    </Grid>
  }, [navEntryId, setShowJson, sectionIsEditable, parentRelation, lane, history, handleArchiveChanged, section])

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

  const subSectionsToRender = sub_sections
    .filter(subSectionDef => section[subSectionDef.name] || config.showAllDefined || sectionIsEditable)
    .filter(filter)
  const subSectionCompartment = (
    <Compartment title="sub sections">
      {subSectionsToRender
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
      {def.m_annotations?.plot && <SectionPlots sectionDef={def} section={section}/>}
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
      {def.m_annotations?.plot && <SectionPlots sectionDef={def} section={section}/>}
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
    const sectionDef = subSectionDef.sub_section
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
  const values = useMemo(() => lane.adaptor.obj[quantityDef.name].map(() => null), [lane.adaptor.obj, quantityDef.name])
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

const usePlotStyle = makeStyles(theme => ({
  error: {
    margin: theme.spacing(1),
    minWidth: 300
  }
}))

const XYPlot = React.memo(function XYPlot({plot, section, sectionDef, title}) {
  const classes = usePlotStyle()
  const theme = useTheme()
  const units = useUnits()
  const xAxis = plot.x || plot['x_axis'] || plot['xAxis']
  const yAxis = plot.y || plot['y_axis'] || plot['yAxis']

  const [data, layout] = useMemo(() => {
    const Y = Array.isArray(yAxis) ? yAxis : [yAxis]
    const nLines = Y.length
    const toUnit = path => {
      const relativePath = path.replace('./', '')
      const resolvedQuantityDef = resolveRef(relativePath, sectionDef)
      if (resolvedQuantityDef === undefined || resolvedQuantityDef === null) {
        throw new Error(`Could not resolve the path ${path}`)
      }
      const value = resolveRef(relativePath, section)
      if (value === undefined || value === null) {
        throw new Error(`Could not resolve the data for ${path}`)
      }
      const unit = resolvedQuantityDef?.unit
      if (unit) {
        const quantity = new Q(value, unit).toSystem(units)
        return [quantity.value(), quantity.label()]
      } else {
        return [value, unit]
      }
    }

    let xValues, xUnit
    try {
      [xValues, xUnit] = toUnit(xAxis)
    } catch (e) {
      return [{error: e.message}, undefined]
    }
    const xPath = xAxis.split('/')
    const xLabel = titleCase(xPath[xPath.length - 1])

    const lines = getLineStyles(nLines, theme).map(line => {
      return {type: 'scatter',
        mode: 'lines',
        line: line}
    })
    if (plot.lines) {
      Y.forEach((y, index) => {
        _.merge(lines[index], plot.lines[index])
      })
    }

    let data = []
    const yUnits = []
    const yLabels = []
    Y.forEach((y, index) => {
      let yValues, yUnit
      try {
        [yValues, yUnit] = toUnit(y)
      } catch (e) {
        data = {error: e.message}
        return
      }
      const yPath = y.split('/')
      const yLabel = titleCase(yPath[yPath.length - 1])
      const line = {
        name: yLabel,
        x: xValues,
        y: yValues,
        ...lines[index]
      }
      data.push(line)
      yUnits.push(yUnit)
      yLabels.push(yLabel)
    })

    const getColor = index => {
      const line = lines[index]
      if ('mode' in line) {
        if (line.mode === 'lines') {
          return {color: line.line?.color}
        } else if (line.mode === 'markers') {
          return {color: line.marker?.color}
        }
      }
      return {color: '#000000'}
    }

    const sameUnit = yUnits.every(unit => unit === yUnits[0])

    const layout = {
      xaxis: {
        title: xUnit ? `${xLabel} (${xUnit})` : xLabel
      },
      yaxis: {
        title: sameUnit ? (yUnits[0] ? `${titleCase(title)} (${yUnits[0]})` : titleCase(title)) : (yUnits[0] ? `${yLabels[0]} (${yUnits[0]})` : yLabels[0]),
        titlefont: !sameUnit && nLines > 1 ? getColor(0) : undefined,
        tickfont: !sameUnit && nLines > 1 ? getColor(0) : undefined
      },
      showlegend: sameUnit && nLines > 1,
      legend: {
        x: 1,
        y: 1,
        xanchor: 'right'
      }
    }

    if (!sameUnit) {
      Y.forEach((y, index) => {
        const color = getColor(index)
        if (index > 0) {
          layout[`yaxis${index + 1}`] = {
            title: yUnits[index] ? `${yLabels[index]} (${yUnits[index]})` : yLabels[index],
            anchor: 'x',
            overlaying: 'y',
            side: index % 2 === 0 ? 'left' : 'right',
            titlefont: nLines > 1 ? color : undefined,
            tickfont: nLines > 1 ? color : undefined
          }
          data[index]['yaxis'] = `y${index + 1}`
        }
      })
    }

    if (plot.layout) {
      _.merge(layout, plot.layout)
    }

    return [data, layout]
  }, [plot.layout, plot.lines, xAxis, yAxis, section, sectionDef, theme, title, units])

  if ('error' in data) {
    return <Alert
      severity="error"
      className={classes.error}
    >
      {`Error when plotting ${titleCase(plot.label)}: ${data.error}`}
    </Alert>
  }

  return <Box minWidth={500}>
    <Plot
      data={data}
      layout={layout}
      floatTitle={title}
      fixedMargins={true}
      config={plot.config}
    />
  </Box>
})
XYPlot.propTypes = {
  plot: PropTypes.object.isRequired,
  sectionDef: PropTypes.object.isRequired,
  section: PropTypes.object,
  title: PropTypes.string
}

export const SectionPlots = React.memo(function SectionPlots({section, sectionDef}) {
  const classes = usePlotStyle()
  const plot = sectionDef.m_annotations?.plot
  const [selected, setSelected] = useState([0])
  const plots = useMemo(() => {
    const plots = (Array.isArray(plot) ? [...plot] : [{...plot}])
    plots.forEach(plot => {
      if (!('label' in plot)) {
        const yAxis = plot.y || plot['y_axis'] || plot['yAxis']
        const pathParts = Array.isArray(yAxis) ? ['unnamed'] : yAxis.split('/')
        plot.label = pathParts[pathParts.length - 1]
      }
    })
    return plots
  }, [plot])

  useEffect(() => {
    setSelected([0])
  }, [plots.length])

  if (plots.length < 1 || selected.find(index => index >= plots.length)) {
    return ''
  }

  const errorMessage = `Unexpected error when plotting ${titleCase(plots[0].label || '')}`

  return <Compartment title="plot">
    <Box minWidth={500}>
      {plots.length > 1 && <TextField
        select variant='filled' size='small' fullWidth label={'Shown plots'} style={{marginBottom: 2}}
        SelectProps={{
          multiple: true,
          value: selected,
          onChange: event => setSelected(event.target.value),
          renderValue: newSelected => newSelected.map(index => titleCase(plots[index].label))?.join(', ')
        }}
      >
        {plots.map((plot, index) => (
          <MenuItem key={index} value={index}>
            <Checkbox
              checked={selected.findIndex(selectedIndex => selectedIndex === index) >= 0}
            />
            {titleCase(plot.label)}
          </MenuItem>
        ))}
      </TextField>}
      {selected.map(index => plots?.[index])
        ?.map((plot, index) => (
          <ErrorHandler
            key={index}
            className={classes.error}
            message={errorMessage}
          >
            <XYPlot
              sectionDef={sectionDef}
              section={section}
              plot={plot}
              title={plot.label}
            />
          </ErrorHandler>
        ))
      }
    </Box>
  </Compartment>
})
SectionPlots.propTypes = {
  sectionDef: PropTypes.object.isRequired,
  section: PropTypes.object
}

function Quantity({value, def}) {
  const {prev} = useLane()
  return <Content>
    <ArchiveTitle def={def} data={value} kindLabel="value" />
    {def.m_annotations?.plot && (
      <Compartment title="plot">
        <XYPlot
          section={prev.adaptor.obj}
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
  const {update, adaptor: {context}} = useContext(laneContext)
  const upload_id = context.upload_id
  context.resources = context.resources || {}
  const resources = context.resources

  useEffect(() => {
    const url = value.split('#')[0]
    if (resources?.[url]) {
      setLoading(false)
      return
    }

    if (!(url.startsWith('../upload/archive/') && upload_id)) {
      setLoading(false)
      return
    }

    api.get(`uploads/${upload_id}/${url.slice('../upload/'.length)}`)
      .then(response => {
        resources[url] = response.data.archive
        update()
      })
      .catch(raiseError)
  }, [api, upload_id, resources, raiseError, setLoading, update, value])

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
