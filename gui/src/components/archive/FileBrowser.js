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
import React, { useContext, useMemo, useCallback, useState, useRef, useEffect } from 'react'
import { useHistory } from 'react-router-dom'
import PropTypes from 'prop-types'
import { makeStyles, Typography, IconButton, Box, Grid, Button, Tooltip, TextField,
  Dialog, DialogContent, DialogContentText } from '@material-ui/core'
import DialogActions from '@material-ui/core/DialogActions'
import { useDataStore, useEntryStoreObj } from '../DataStore'
import Browser, { Item, Content, Adaptor, browserContext, laneContext, Title, Compartment } from './Browser'
import { useApi } from '../api'
import UploadIcon from '@material-ui/icons/CloudUpload'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import FileCopyIcon from '@material-ui/icons/FileCopy'
import CreateNewFolderIcon from '@material-ui/icons/CreateNewFolder'
import DeleteIcon from '@material-ui/icons/Delete'
import NavigateIcon from '@material-ui/icons/ArrowForward'
import FolderIcon from '@material-ui/icons/FolderOutlined'
import FileIcon from '@material-ui/icons/InsertDriveFileOutlined'
import ReloadIcon from '@material-ui/icons/Replay'
import RecognizedFileIcon from '@material-ui/icons/InsertChartOutlinedTwoTone'
import { useDropzone } from 'react-dropzone'
import { fromEvent } from 'file-selector'
import Download from '../entry/Download'
import Quantity from '../Quantity'
import FilePreview from './FilePreview'
import { useEntryStore } from '../entry/EntryContext'
import { archiveAdaptorFactory } from './ArchiveBrowser'
import NorthLaunchButton from '../north/NorthLaunchButton'
import { useTools } from '../north/NorthPage'
import { EntryButton } from '../nav/Routes'
import { useErrors } from '../errors'
import { apiBase } from '../../config'
import { parseNomadUrl, refType, urlJoin, urlEncodePath, systemMetainfoUrl, createEntryUrl } from '../../utils'

const FileBrowser = React.memo(({uploadUrl, rootTitle, highlightedItem = null}) => {
  const adaptor = useMemo(() => {
    return new RawDirectoryAdaptor(uploadUrl, rootTitle, highlightedItem)
  }, [uploadUrl, rootTitle, highlightedItem])

  return <Browser adaptor={adaptor} />
})
FileBrowser.propTypes = {
  uploadUrl: PropTypes.string.isRequired,
  rootTitle: PropTypes.string.isRequired,
  highlightedItem: PropTypes.string
}
export default FileBrowser

class RawDirectoryAdaptor extends Adaptor {
  constructor(uploadUrl, title, highlightedItem) {
    super()
    const parsedUrl = parseNomadUrl(uploadUrl)
    if (parsedUrl.type !== refType.upload) throw new Error(`Expected an upload url, got ${uploadUrl}`)
    if (!parsedUrl.isResolved) throw new Error(`Absolute url required, got ${uploadUrl}`)
    this.deploymentUrl = parsedUrl.deploymentUrl
    this.uploadUrl = uploadUrl
    this.uploadId = parsedUrl.uploadId
    this.path = parsedUrl.path
    this.title = title
    this.highlightedItem = highlightedItem
    this.editable = undefined
    this.data = undefined
    this.timestamp = undefined
    this.initialized = false
    this.dependencies = new Set([this.uploadId])
    this.api = undefined
    this.page_size = 40
    this.lastPage = undefined
  }
  depends() {
    return this.dependencies
  }
  async initialize(api, dataStore) {
    this.api = api
    this.lastPage = 1
    const uploadStoreObj = await dataStore.getUploadAsync(this.deploymentUrl, this.uploadId, true, false)
    this.timestamp = uploadStoreObj.upload?.complete_time
    this.editable = uploadStoreObj.isEditable
    await this.fetchData()
  }
  cleanup() {
    delete this.data
  }
  async fetchData() {
    const encodedPath = urlEncodePath(this.path)
    if (this.deploymentUrl !== apiBase) throw new Error('Fetching directory data from external source is not yet supported')
    const response = await this.api.get(`/uploads/${this.uploadId}/rawdir/${encodedPath}?include_entry_info=true&page_size=${this.page_size}&page=${this.lastPage}`)
    const elementsByName = this.data?.elementsByName || {}
    response.directory_metadata.content.forEach(element => { elementsByName[element.name] = element })
    const directory_metadata = this.data?.directory_metadata
      ? this?.data?.directory_metadata.concat(response.directory_metadata.content)
      : response.directory_metadata.content
    const total = response.pagination.total
    this.data = {directory_metadata, elementsByName, total}
  }
  async itemAdaptor(key) {
    if (key === '_mainfile') {
      key = this.highlightedItem
    }
    const extendedUrl = urlJoin(this.uploadUrl, encodeURIComponent(key))
    let element = this.data?.elementsByName[key]
    while (!element && this.lastPage * this.page_size < this.data.total) {
      const data = await this.scrollToNextPage()
      element = data?.elementsByName[key]
    }
    if (element) {
      if (element.is_file) {
        return new RawFileAdaptor(extendedUrl, element, this.editable)
      } else {
        return new RawDirectoryAdaptor(extendedUrl, key, null)
      }
    }
    throw new Error('Bad path: ' + key)
  }

  async scrollToNextPage() {
    if (this.data?.total && this.lastPage * this.page_size < this.data.total) {
      this.lastPage = this.lastPage + 1
      await this.fetchData()
      return this.data
    }
  }

  async onScrollToEnd(event, updateLane) {
    const data = await this.scrollToNextPage()
    if (data) {
      updateLane()
    }
  }

  async onRendered(event, updateLane) {
    if (event.target.scrollHeight * 95 / 100 <= event.target.clientHeight) {
      const data = await this.scrollToNextPage()
      if (data) {
        updateLane()
      }
    }
  }

  render() {
    return <RawDirectoryContent
      deploymentUrl={this.deploymentUrl} uploadId={this.uploadId} path={this.path}
      title={this.title} highlightedItem={this.highlightedItem}
      editable={this.editable}/>
  }
}

export const CustomDropZone = React.memo((props) => {
  const {children, onDrop, classNameNormal, classNameActive, ...otherProps} = props
  const [active, setActive] = useState(false)
  const dragTargetsRef = useRef([])

  const onDragEnter = useCallback(e => {
    e.preventDefault()
    e.stopPropagation()
    dragTargetsRef.current.push(e.target)
    setActive(true)
  }, [setActive])

  const onDragLeave = useCallback(e => {
    e.preventDefault()
    e.stopPropagation()
    const targets = dragTargetsRef.current
    const targetIdx = targets.indexOf(e.target)
    if (targetIdx !== -1) {
      targets.splice(targetIdx, 1)
      dragTargetsRef.current = targets
    }
    if (!targets.length) {
      setActive(false)
    }
  }, [setActive])

  const onDragOver = useCallback(e => {
    e.preventDefault()
    e.stopPropagation()
  }, [])

  const onDropped = useCallback(e => {
    e.preventDefault()
    e.stopPropagation()
    setActive(false)
    dragTargetsRef.current = []
    onDrop(e)
  }, [onDrop, setActive])

  return <div
    className={active ? classNameActive : classNameNormal}
    onDragEnter={onDragEnter}
    onDragLeave={onDragLeave}
    onDragOver={onDragOver}
    onDrop={onDropped}
    {...otherProps}
  >
    {children}
  </div>
})
CustomDropZone.propTypes = {
  children: PropTypes.any,
  onDrop: PropTypes.func,
  classNameNormal: PropTypes.string.isRequired,
  classNameActive: PropTypes.string.isRequired
}

const useRawDirectoryContentStyles = makeStyles(theme => ({
  dropzoneLane: {
    width: '100%',
    minHeight: '100%'
  },
  dropzoneTop: {
    width: '100%',
    margin: 0,
    padding: 0
  },
  dropzoneActive: {
    width: '100%',
    minHeight: '100%',
    backgroundColor: theme.palette.grey[300]
  }
}))
function RawDirectoryContent({deploymentUrl, uploadId, path, title, highlightedItem, editable}) {
  const classes = useRawDirectoryContentStyles()
  const dataStore = useDataStore()
  const entryPageMainFile = useEntryStore()?.metadata?.mainfile // Will be set if we're on an entry page
  const browser = useContext(browserContext)
  const lane = useContext(laneContext)
  const history = useHistory()
  const encodedPath = urlEncodePath(path)
  const { api } = useApi()
  const [openConfirmDeleteDirDialog, setOpenConfirmDeleteDirDialog] = useState(false)
  const [openCreateDirDialog, setOpenCreateDirDialog] = useState(false)
  const [openCopyMoveDialog, setOpenCopyMoveDialog] = useState(false)
  const [openDuplicateDialog, setOpenDuplicateDialog] = useState(false)
  const [fileName, setFileName] = useState('')
  const copyFileName = useRef()
  const createDirName = useRef()
  const { raiseError } = useErrors()

  const refreshIfNeeded = useCallback(async (oldStoreObj, newStoreObj) => {
    // Called when a store update occurs. Note, adaptors can be cached while the browser & lanes
    // components are unmounted.
    if (newStoreObj.hasUpload) {
      const oldTimestamp = lane.adaptor.timestamp // Timestamp from the store the last time we called initialize
      const newTimestamp = newStoreObj.upload?.complete_time // Current timestamp from the store
      if (!lane.adaptor.initialized) {
        // No need to initiate a new refresh, just update the adaptor timestamp
        lane.adaptor.timestamp = newTimestamp
        lane.adaptor.initialized = true
      } else if (newTimestamp !== oldTimestamp && !newStoreObj.isProcessing) {
        // Reprocessed
        browser.invalidateLanesWithDependency(uploadId)
      }
    }
  }, [browser, lane, uploadId])

  useEffect(() => {
    refreshIfNeeded(undefined, dataStore.getUpload(deploymentUrl, uploadId))
    return dataStore.subscribeToUpload(deploymentUrl, uploadId, refreshIfNeeded, true, false)
  }, [dataStore, deploymentUrl, uploadId, refreshIfNeeded])

  const handleDropFiles = useCallback((files) => {
    // Handles dropping files (not links)
    const formData = new FormData() // eslint-disable-line no-undef
    for (const file of files) {
      formData.append('file', file)
    }
    api.put(`/uploads/${uploadId}/raw/${encodedPath}`, formData)
      .then(response => dataStore.updateUpload(deploymentUrl, uploadId, {upload: response.data}))
      .catch(error => raiseError(error))
  }, [deploymentUrl, uploadId, encodedPath, raiseError, api, dataStore])

  const handleDrop = useCallback(async (e) => {
    // Handles dropping files and links
    const files = await fromEvent(e)
    const _filePath = e.dataTransfer?.getData('URL')
    if (files.length) { // files are being transferred
      handleDropFiles(files)
    } else if (_filePath) {
      if (_filePath.includes('/files/') &&
      _filePath.includes(history.location.pathname.split('files')[0])) {
        setFileName(_filePath.slice(_filePath.indexOf('files')).split('/').slice(1).join('/'))
        setOpenCopyMoveDialog(true)
      }
    }
  }, [history.location.pathname, handleDropFiles, setFileName, setOpenCopyMoveDialog])

  const handleCopyMoveFile = (e) => {
    setOpenCopyMoveDialog(false)
    setOpenDuplicateDialog(false)
    api.put(`/uploads/${uploadId}/raw/${encodedPath}?copy_or_move=${e.moveFile}&copy_or_move_source_path=${fileName}&file_name=${copyFileName.current.value}`)
      .then(response => dataStore.updateUpload(deploymentUrl, uploadId, {upload: response.data}))
      .catch(error => raiseError(error))
  }

  const handleCreateDir = () => {
    setOpenCreateDirDialog(false)
    const dirName = createDirName.current.value
    if (dirName) {
      const fullPath = urlJoin(encodedPath, encodeURIComponent(dirName))
      api.post(`/uploads/${uploadId}/raw-create-dir/${fullPath}`)
        .then(response => dataStore.updateUpload(deploymentUrl, uploadId, {upload: response.data}))
        .catch(raiseError)
    }
  }

  const handleDeleteDir = () => {
    setOpenConfirmDeleteDirDialog(false)
    api.delete(`/uploads/${uploadId}/raw/${encodedPath}`)
      .then(response => {
        if (typeof mainfile === 'string' && (path === '' || entryPageMainFile === path || entryPageMainFile.startsWith(path + '/'))) {
          // This will delete the current entry - go to upload overview page
          history.push(`/user/uploads/upload/id/${uploadId}`)
        } else {
          const gotoLane = lane.index > 0 ? lane.index - 1 : 0
          history.push(browser.lanes.current[gotoLane].path)
        }
        dataStore.updateUpload(deploymentUrl, uploadId, {upload: response.data})
      })
      .catch(raiseError)
  }

  const { getRootProps, getInputProps } = useDropzone({onDrop: handleDropFiles})

  if (!lane.adaptor.data) {
    return <Content key={path}><Typography>loading ...</Typography></Content>
  } else {
    // Data loaded
    const downloadUrl = `uploads/${uploadId}/raw/${encodedPath}?compress=true` // TODO: deploymentUrl need to be considered for external uploads
    return (
      <CustomDropZone onDrop={handleDrop} classNameNormal={classes.dropzoneLane} classNameActive={classes.dropzoneActive}>
        <Content key={path}>
          <Title
            title={title}
            label="folder"
            tooltip={path}
            actions={
              <Grid container justifyContent="space-between" wrap="nowrap" spacing={1}>
                <Grid item>
                  <IconButton size="small" onClick={() => setOpenDuplicateDialog(true)} disabled={!editable || fileName === ""}>
                    <Tooltip title="duplicate file">
                      <FileCopyIcon />
                    </Tooltip>
                  </IconButton>
                  <Dialog
                      open={openDuplicateDialog}
                      onClose={() => setOpenDuplicateDialog(false)}
                    >
                      <DialogContent>
                        <DialogContentText>Duplicate <b>{fileName}</b>:</DialogContentText>
                        <TextField
                          fullWidth
                          placeholder='Provide a name'
                          autoFocus
                          inputRef={copyFileName}
                          defaultValue={`(Copy) ${fileName.split('/').splice(-1)}`}
                        />
                      </DialogContent>
                      <DialogActions>
                        <Button onClick={() => setOpenDuplicateDialog(false)}>Cancel</Button>
                        <Button onClick={(e) => handleCopyMoveFile({...e, moveFile: 'copy'})}>Duplicate</Button>
                      </DialogActions>
                    </Dialog>
                </Grid>
                <Grid item>
                  <IconButton size="small" onClick={() => browser.invalidateLanesFromIndex(lane.index)}>
                    <Tooltip title="reload directory contents">
                      <ReloadIcon/>
                    </Tooltip>
                  </IconButton>
                </Grid>
                <Grid item>
                  <Download
                    component={IconButton} disabled={false} size="small"
                    tooltip="download this folder"
                    url={downloadUrl}
                  >
                    <DownloadIcon />
                  </Download>
                </Grid>
                {
                  editable &&
                    <Grid item>
                      <div {...getRootProps()}>
                        <input {...getInputProps()} />
                        <IconButton size="small">
                          <Tooltip title="upload to this folder (click or drop files)">
                            <UploadIcon/>
                          </Tooltip>
                        </IconButton>
                      </div>
                    </Grid>
                }
                {
                  editable &&
                    <Grid item>
                      <IconButton size="small" onClick={() => setOpenCreateDirDialog(true)}>
                        <Tooltip title="create new folder">
                          <CreateNewFolderIcon />
                        </Tooltip>
                      </IconButton>
                      <Dialog
                        open={openCreateDirDialog}
                        onClose={() => setOpenCreateDirDialog(false)}
                      >
                        <DialogContent>
                          <DialogContentText>Name of new directory:</DialogContentText>
                          <TextField
                            autoFocus
                            inputRef={createDirName}
                            onKeyDown={(event) => {
                              if (event.key === 'Enter') {
                                handleCreateDir()
                              }
                            }}/>
                        </DialogContent>
                        <DialogActions>
                          <Button onClick={() => setOpenCreateDirDialog(false)}>Cancel</Button>
                          <Button onClick={() => handleCreateDir()}>OK</Button>
                        </DialogActions>
                      </Dialog>
                    </Grid>
                }
                {
                  editable &&
                    <Dialog
                      open={openCopyMoveDialog}
                      onClose={() => setOpenCopyMoveDialog(false)}
                    >
                      <DialogContent>
                        <DialogContentText>Copy/move <b>{fileName}</b>:</DialogContentText>
                        <TextField
                          fullWidth
                          placeholder='Provide a name'
                          autoFocus
                          inputRef={copyFileName}
                          defaultValue={fileName.split('/').splice(-1)}
                        />
                      </DialogContent>
                      <DialogActions>
                        <Button onClick={() => setOpenCopyMoveDialog(false)}>Cancel</Button>
                        <Button onClick={(e) => handleCopyMoveFile({...e, moveFile: 'copy'})}>Copy</Button>
                        <Button onClick={(e) => handleCopyMoveFile({...e, moveFile: 'move'})}>Move</Button>
                      </DialogActions>
                    </Dialog>
                }
                {
                  editable &&
                    <Grid item>
                      <IconButton size="small" onClick={() => setOpenConfirmDeleteDirDialog(true)}>
                        <Tooltip title="delete this folder">
                          <DeleteIcon />
                        </Tooltip>
                      </IconButton>
                      <Dialog
                        open={openConfirmDeleteDirDialog}
                        onClose={() => setOpenConfirmDeleteDirDialog(false)}
                      >
                        <DialogContent>
                          <DialogContentText>Really delete the directory <b>{path}</b>?</DialogContentText>
                        </DialogContent>
                        <DialogActions>
                          <Button onClick={() => setOpenConfirmDeleteDirDialog(false)} autoFocus>Cancel</Button>
                          <Button onClick={handleDeleteDir}>OK</Button>
                        </DialogActions>
                      </Dialog>
                    </Grid>
                }
              </Grid>
            }
          />
          <Compartment>
            {
              lane.adaptor.data.directory_metadata.map(element => (
                <span key={element.name} onClick={() => setFileName(element.name)}>
                  <Item
                    icon={element.is_file ? (element.parser_name ? RecognizedFileIcon : FileIcon) : FolderIcon}
                    itemKey={element.name} key={urlJoin(path, element.name)}
                    highlighted={element.name === highlightedItem}
                    chip={element.parser_name && element.parser_name.replace('parsers/', '').replace('archive', 'nomad')}
                  >
                    {element.name}
                  </Item>
                </span>
              ))
            }
          </Compartment>
        </Content>
      </CustomDropZone>
    )
  }
}
RawDirectoryContent.propTypes = {
  deploymentUrl: PropTypes.string.isRequired,
  uploadId: PropTypes.string.isRequired,
  path: PropTypes.string.isRequired,
  title: PropTypes.string.isRequired,
  highlightedItem: PropTypes.string,
  editable: PropTypes.bool.isRequired
}

export class RawFileAdaptor extends Adaptor {
  /**
   * Constructs an adaptor for a file lane. Note, data is optional, if not provided it will
   * be fetched in initialize.
   */
  constructor(uploadUrl, data, editable) {
    super()
    const parsedUrl = parseNomadUrl(uploadUrl)
    if (parsedUrl.type !== refType.upload) throw new Error(`Expected an upload url, got ${uploadUrl}`)
    if (!parsedUrl.isResolved) throw new Error(`Absolute url required, got ${uploadUrl}`)
    this.deploymentUrl = parsedUrl.deploymentUrl
    this.uploadId = parsedUrl.uploadId
    this.path = parsedUrl.path
    this.data = data
    this.editable = editable
    this.dependencies = new Set([this.uploadId])
  }
  depends() {
    return this.dependencies
  }
  async initialize(api, dataStore) {
    this.dataStore = dataStore
    if (!this.data) {
      const response = await api.get(`uploads/${this.uploadId}/rawdir/${this.path}`)
      this.data = response.file_metadata
    }
  }
  async itemAdaptor(key) {
    if (key === 'archive') {
      const archiveUrl = createEntryUrl(this.deploymentUrl, this.uploadId, this.data.entry_id)
      const {archive} = await this.dataStore.getEntryAsync(this.deploymentUrl, this.data.entry_id, false, '*')
      const metainfo = await this.dataStore.getMetainfoAsync(systemMetainfoUrl)
      const rootSectionDef = metainfo.getEntryArchiveDefinition()
      return archiveAdaptorFactory(archiveUrl, archive, rootSectionDef)
    } else if (key === 'preview') {
      return new FilePreviewAdaptor(this.uploadId, this.path, this.data)
    }
  }
  render() {
    return <RawFileContent
      deploymentUrl={this.deploymentUrl} uploadId={this.uploadId} path={this.path}
      data={this.data} editable={this.editable} key={this.path}/>
  }
}

class FilePreviewAdaptor extends Adaptor {
  constructor(uploadId, path, data) {
    super()
    this.uploadId = uploadId
    this.path = path
    this.data = data
  }
  render() {
    return <Box minWidth={'300px'} maxWidth={'1200px'} width="100%">
      <FilePreview uploadId={this.uploadId} path={this.path} size={this.data?.size}/>
    </Box>
  }
}

function RawFileContent({deploymentUrl, uploadId, path, data, editable}) {
  const entryPageMainFile = useEntryStore()?.metadata?.mainfile // Will be set if we're on an entry page
  const browser = useContext(browserContext)
  const lane = useContext(laneContext)
  const history = useHistory()
  const dataStore = useDataStore()
  const { api } = useApi()
  const { raiseError } = useErrors()
  const [openConfirmDeleteFileDialog, setOpenConfirmDeleteFileDialog] = useState(false)
  const encodedPath = urlEncodePath(path)
  const downloadUrl = `uploads/${uploadId}/raw/${encodedPath}?ignore_mime_type=true` // TODO: deploymentUrl need to be considered for external uploads
  const allNorthTools = useTools()
  const applicableNorthTools = useMemo(() => {
    if (!allNorthTools || typeof allNorthTools !== 'object' || Object.keys(allNorthTools).length < 1) {
      return []
    }
    const fileExtension = path.split('.').pop().toLowerCase()
    return Object.keys(allNorthTools)
      .filter(key => {
        const tool = allNorthTools[key]
        return tool?.file_extensions && tool.file_extensions.includes(fileExtension)
      })
  }, [allNorthTools, path])

  // If we are stepping into an archive, subscribe to the store for updates on this archive
  const archiveSelected = data.entry_id && lane.next?.key === 'archive'
  const entryStoreObj = useEntryStoreObj(archiveSelected ? deploymentUrl : null, data.entry_id, false, '*')

  // Invalidate subsequent lanes if the archive is updated
  useEffect(() => {
    if (entryStoreObj?.archive && lane.next?.adaptor?.obj) {
      if (entryStoreObj.archive !== lane.next.adaptor?.obj) {
        browser.invalidateLanesFromIndex(lane.index)
      }
    }
  }, [browser, lane, entryStoreObj])

  const handleDeleteFile = () => {
    setOpenConfirmDeleteFileDialog(false)
    api.delete(`/uploads/${uploadId}/raw/${encodedPath}`)
      .then(response => {
        if (path === entryPageMainFile) {
          // Deleting the main entry - go to upload overview page
          history.push(`/user/uploads/upload/id/${uploadId}`)
        } else {
          const gotoLane = lane.index > 0 ? lane.index - 1 : 0
          history.push(browser.lanes.current[gotoLane].path)
        }
        dataStore.updateUpload(deploymentUrl, uploadId, {upload: response.data})
      })
      .catch(raiseError)
  }

  // A nicer, human-readable size string
  let niceSize, unit, factor
  if (data.size > 1e9) {
    [unit, factor] = ['GB', 1e9]
  } else if (data.size > 1e6) {
    [unit, factor] = ['MB', 1e6]
  } else if (data.size > 1e3) {
    [unit, factor] = ['kB', 1e3]
  }
  if (unit) {
    if (data.size / factor > 100) {
      // No decimals
      niceSize = `${Math.round(data.size / factor)} ${unit} (${data.size} bytes)`
    } else {
      // One decimal
      niceSize = `${Math.round(data.size * 10 / factor) / 10} ${unit} (${data.size} bytes)`
    }
  } else {
    // Unit = bytes
    niceSize = `${data.size} bytes`
  }

  return (
    <Content
      key={path}
      display="flex" flexDirection="column" height="100%"
      paddingBottom={0} maxWidth="initial"
    >
      <Title
        title={data.name}
        label="file"
        tooltip={path}
        actions={
          <Grid container justifyContent="space-between" wrap="nowrap" spacing={1}>
            {data.entry_id && (
              <Grid item>
                <EntryButton entryId={data.entry_id} component={IconButton} size="small">
                  <NavigateIcon />
                </EntryButton>
              </Grid>
            )}
            <Grid item>
              <Download
                component={IconButton} disabled={false} size="small"
                tooltip="download this file"
                url={downloadUrl}
              >
                <DownloadIcon />
              </Download>
            </Grid>
            {
              editable &&
                <Grid item>
                  <IconButton size="small" onClick={() => setOpenConfirmDeleteFileDialog(true)}>
                    <Tooltip title="delete this file">
                      <DeleteIcon />
                    </Tooltip>
                  </IconButton>
                  <Dialog
                    open={openConfirmDeleteFileDialog}
                    onClose={() => setOpenConfirmDeleteFileDialog(false)}
                  >
                    <DialogContent>
                      <DialogContentText>Really delete the file <b>{data.name}</b>?</DialogContentText>
                    </DialogContent>
                    <DialogActions>
                      <Button onClick={() => setOpenConfirmDeleteFileDialog(false)} autoFocus>Cancel</Button>
                      <Button onClick={handleDeleteFile}>OK</Button>
                    </DialogActions>
                  </Dialog>
                </Grid>
            }
          </Grid>
        }
      />
      <Compartment>
        <Quantity quantity="filesize" data={{filesize: niceSize}} label="file size" />
        {data.parser_name && (
          <Quantity quantity="parser" data={{parser: data.parser_name}} />
        )}
        {data.entry_id && (
          <Quantity quantity="entryId" data={{entryId: data.entry_id}} noWrap withClipboard />
        )}
      </Compartment>
      <Compartment>
        {data.entry_id && (
          <Item itemKey="archive">
            <Typography>processed data</Typography>
          </Item>
        )}
        <Item itemKey="preview">
          <Typography>preview</Typography>
        </Item>
        {applicableNorthTools.length > 0 && (
          <NorthLaunchButton uploadId={uploadId} path={path} tools={applicableNorthTools} />
        )}
      </Compartment>
    </Content>)
}
RawFileContent.propTypes = {
  deploymentUrl: PropTypes.string.isRequired,
  uploadId: PropTypes.string.isRequired,
  path: PropTypes.string.isRequired,
  data: PropTypes.object.isRequired,
  editable: PropTypes.bool.isRequired
}
