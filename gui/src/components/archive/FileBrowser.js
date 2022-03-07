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
import React, { useContext, useMemo, useState } from 'react'
import { useHistory } from 'react-router-dom'
import PropTypes from 'prop-types'
import { makeStyles, Typography, IconButton, Box, Grid, Button, Tooltip,
  Dialog, DialogContent, DialogContentText } from '@material-ui/core'
import DialogActions from '@material-ui/core/DialogActions'
import Browser, { Item, Content, Adaptor, browserContext, laneContext, Title, Compartment } from './Browser'
import { useApi } from '../api'
import UploadIcon from '@material-ui/icons/CloudUpload'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import DeleteIcon from '@material-ui/icons/Delete'
import NavigateIcon from '@material-ui/icons/MoreHoriz'
import FolderIcon from '@material-ui/icons/FolderOutlined'
import FileIcon from '@material-ui/icons/InsertDriveFileOutlined'
import ReloadIcon from '@material-ui/icons/Replay'
import RecognizedFileIcon from '@material-ui/icons/InsertChartOutlinedTwoTone'
import Dropzone from 'react-dropzone'
import Download from '../entry/Download'
import Quantity from '../Quantity'
import FilePreview from './FilePreview'
import { archiveAdaptorFactory } from './ArchiveBrowser'
import H5Web from '../visualization/H5Web'
import NorthLaunchButton from '../north/NorthLaunchButton'
import { useTools } from '../north/NorthPage'
import { EntryButton } from '../nav/Routes'

const FileBrowser = React.memo(({uploadId, path, rootTitle, highlightedItem = null, editable = false}) => {
  const adaptor = new RawDirectoryAdaptor(uploadId, path, rootTitle, highlightedItem, editable)
  return <Browser adaptor={adaptor} />
})
FileBrowser.propTypes = {
  uploadId: PropTypes.string.isRequired,
  path: PropTypes.string.isRequired,
  rootTitle: PropTypes.string.isRequired,
  highlightedItem: PropTypes.string,
  editable: PropTypes.bool
}
export default FileBrowser

class RawDirectoryAdaptor extends Adaptor {
  constructor(uploadId, path, title, highlightedItem, editable = false) {
    super()
    this.uploadId = uploadId
    this.path = path
    this.title = title
    this.highlightedItem = highlightedItem
    this.editable = editable
    this.data = undefined // Will be set by RawDirectoryContent component when loaded
  }
  needToFetchData() {
    return this.data === undefined
  }
  async fetchData(api) {
    if (this.data === undefined) {
      const encodedPath = this.path.split('/').map(segment => encodeURIComponent(segment)).join('/')
      const response = await api.get(`/uploads/${this.uploadId}/rawdir/${encodedPath}?include_entry_info=true&page_size=500`)
      const elementsByName = {}
      response.directory_metadata.content.forEach(element => { elementsByName[element.name] = element })
      this.data = {response, elementsByName}
    }
  }
  onFilesUpdated(uploadId, path) {
    if (uploadId === this.uploadId && this.path.startsWith(path)) {
      this.data = undefined
    }
  }
  itemAdaptor(key) {
    const ext_path = this.path ? this.path + '/' + key : key
    const element = this.data.elementsByName[key]
    if (element) {
      if (element.is_file) {
        return new RawFileAdaptor(this.uploadId, ext_path, element, this.editable)
      } else {
        return new RawDirectoryAdaptor(this.uploadId, ext_path, key, null, this.editable)
      }
    }
    throw new Error('Bad path: ' + key)
  }
  render() {
    return <RawDirectoryContent
      uploadId={this.uploadId} path={this.path} title={this.title} highlightedItem={this.highlightedItem}
      editable={this.editable}/>
  }
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
    backgroundColor: theme.palette.grey[300]
  }
}))
function RawDirectoryContent({uploadId, path, title, highlightedItem, editable}) {
  const classes = useRawDirectoryContentStyles()
  const browser = useContext(browserContext)
  const lane = useContext(laneContext)
  const history = useHistory()
  const encodedPath = path.split('/').map(segment => encodeURIComponent(segment)).join('/')
  const { api } = useApi()
  const [openConfirmDeleteDirDialog, setOpenConfirmDeleteDirDialog] = useState(false)

  const handleDrop = (files) => {
    if (!files[0]?.name) {
      return // Not dropping a file, but something else. Ignore.
    }
    const formData = new FormData() // eslint-disable-line no-undef
    formData.append('file', files[0])
    browser.blockUntilProcessed({
      uploadId: uploadId,
      apiCall: api.put(`/uploads/${uploadId}/raw/${encodedPath}`, formData, {
        onUploadProgress: (progressEvent) => {
          // TODO: would be nice to show progress somehow
        }
      }),
      apiCallText: 'Uploading file',
      onSuccess: () => {
        browser.lanes.current.forEach(lane => lane.adaptor?.onFilesUpdated(uploadId, path))
        browser.update()
      }
    })
  }

  const handleDeleteDir = () => {
    setOpenConfirmDeleteDirDialog(false)
    browser.blockUntilProcessed({
      uploadId: uploadId,
      apiCall: api.delete(`/uploads/${uploadId}/raw/${encodedPath}`),
      apiCallText: 'Deleting directory',
      onSuccess: () => {
        const segments = path.split('/')
        const parentPath = segments.slice(0, segments.length - 1).join('/')
        browser.lanes.current.forEach(lane => lane.adaptor?.onFilesUpdated(uploadId, parentPath))
        if (lane.index > 0) {
          history.push(browser.lanes.current[lane.index - 1].path)
        }
      }})
  }

  if (lane.adaptor.data === undefined) {
    return <Content key={path}><Typography>loading ...</Typography></Content>
  } else {
    // Data loaded
    const downloadUrl = `uploads/${uploadId}/raw/${encodedPath}?compress=true`
    const segments = path.split('/')
    const lastSegment = segments[segments.length - 1]
    const downloadFilename = `${uploadId}${lastSegment ? ' - ' + lastSegment : ''}.zip`
    return (
      <Dropzone
        disabled={!editable}
        className={classes.dropzoneLane}
        activeClassName={classes.dropzoneActive}
        onDrop={handleDrop} disableClick
      >
        <Content key={path}>
          <Title
            title={title}
            label="folder"
            tooltip={path}
            actions={
              <Grid container justifyContent="space-between" wrap="nowrap" spacing={1}>
                <Grid item>
                  <IconButton size="small" onClick={() => browser.update(lane)}>
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
                    fileName={downloadFilename}
                  >
                    <DownloadIcon />
                  </Download>
                </Grid>
                {
                  editable &&
                    <Grid item>
                      <Dropzone
                        className={classes.dropzoneTop} activeClassName={classes.dropzoneActive}
                        onDrop={handleDrop}
                      >
                        <IconButton size="small">
                          <Tooltip title="upload to this folder (click or drop files)">
                            <UploadIcon/>
                          </Tooltip>
                        </IconButton>
                      </Dropzone>
                    </Grid>
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
              lane.adaptor.data.response.directory_metadata.content.map(element => (
                <Item
                  icon={element.is_file ? (element.parser_name ? RecognizedFileIcon : FileIcon) : FolderIcon}
                  itemKey={element.name} key={path ? path + '/' + element.name : element.name}
                  highlighted={element.name === highlightedItem}
                  chip={element.parser_name && element.parser_name.replace('parsers/', '')}
                >
                  {element.name}
                </Item>
              ))
            }
            {
              lane.adaptor.data.response.pagination.total > 500 &&
                <Typography color="error">Only showing the first 500 rows</Typography>
            }
          </Compartment>
        </Content>
      </Dropzone>)
  }
}
RawDirectoryContent.propTypes = {
  uploadId: PropTypes.string.isRequired,
  path: PropTypes.string.isRequired,
  title: PropTypes.string.isRequired,
  highlightedItem: PropTypes.string,
  editable: PropTypes.bool.isRequired
}

export class RawFileAdaptor extends Adaptor {
  constructor(uploadId, path, data, editable) {
    super()
    this.uploadId = uploadId
    this.path = path
    this.data = data
    this.editable = editable
  }
  async itemAdaptor(key, api) {
    if (key === 'archive') {
      if (!this.data.archive) {
        const response = await api.get(`entries/${this.data.entry_id}/archive`)
        this.data.archive = response.data.archive
      }
      return archiveAdaptorFactory(this.data.archive)
    } else if (key === 'h5web') {
      return new H5WebAdaptor(this.uploadId, this.path)
    }
  }
  render() {
    return <RawFileContent
      uploadId={this.uploadId} path={this.path} data={this.data} editable={this.editable}
      key={this.path}/>
  }
}

class H5WebAdaptor extends Adaptor {
  constructor(uploadId, path) {
    super()
    this.uploadId = uploadId
    this.path = path
  }
  render() {
    return <Box width="80vw" height="100%">
      <H5Web upload_id={this.uploadId} filename={this.path}/>
    </Box>
  }
}

function RawFileContent({uploadId, path, data, editable}) {
  const browser = useContext(browserContext)
  const lane = useContext(laneContext)
  const history = useHistory()
  const { api } = useApi()
  const [openConfirmDeleteFileDialog, setOpenConfirmDeleteFileDialog] = useState(false)
  const encodedPath = path.split('/').map(segment => encodeURIComponent(segment)).join('/')
  const downloadUrl = `uploads/${uploadId}/raw/${encodedPath}?ignore_mime_type=true`
  const allNorthTools = useTools()
  const applicableNorthTools = useMemo(() => {
    const fileExtension = path.split('.').pop().toLowerCase()
    return Object.keys(allNorthTools)
      .filter(key => {
        const tool = allNorthTools[key]
        return tool.file_extensions && tool.file_extensions.includes(fileExtension)
      })
  }, [allNorthTools, path])

  const handleDeleteFile = () => {
    setOpenConfirmDeleteFileDialog(false)
    browser.blockUntilProcessed({
      uploadId: uploadId,
      apiCall: api.delete(`/uploads/${uploadId}/raw/${encodedPath}`),
      apiCallText: 'Deleting file',
      onSuccess: () => {
        const segments = path.split('/')
        const parentPath = segments.slice(0, segments.length - 1).join('/')
        browser.lanes.current.forEach(lane => lane.adaptor?.onFilesUpdated(uploadId, parentPath))
        if (lane.index > 0) {
          history.push(browser.lanes.current[lane.index - 1].path)
        }
      }})
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
      paddingTop={0} paddingBottom={0} maxWidth="initial"
    >
      <Box paddingTop={1}>
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
                  fileName={data.name}
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
      </Box>
      <Compartment>
        <Quantity quantity="filesize" data={{filesize: niceSize}} label="file size" />
        { data.parser_name && <Quantity quantity="parser" data={{parser: data.parser_name}} />}
        { data.parser_name && <Quantity quantity="entryId" data={{entryId: data.entry_id}} noWrap withClipboard />}
      </Compartment>
      {data.entry_id && (
        <Compartment title="more">
          <Item itemKey="archive"><Typography>processed data</Typography></Item>
        </Compartment>
      )}
      {applicableNorthTools.length > 0 && (
        <Compartment title="tools">
          <NorthLaunchButton uploadId={uploadId} path={path} tools={applicableNorthTools} />
        </Compartment>
      )}
      <Compartment title="preview">
        <FilePreview uploadId={uploadId} path={path} size={data.size}/>
      </Compartment>

    </Content>)
}
RawFileContent.propTypes = {
  uploadId: PropTypes.string.isRequired,
  path: PropTypes.string.isRequired,
  data: PropTypes.object.isRequired,
  editable: PropTypes.bool.isRequired
}
