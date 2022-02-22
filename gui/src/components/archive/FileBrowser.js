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
import React, { useContext } from 'react'
import PropTypes from 'prop-types'
import { Typography, IconButton, Box } from '@material-ui/core'
import Browser, { Item, Content, Adaptor, laneContext, Title, Compartment } from './Browser'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import FolderIcon from '@material-ui/icons/FolderOutlined'
import FileIcon from '@material-ui/icons/InsertDriveFileOutlined'
import RecognizedFileIcon from '@material-ui/icons/InsertChartOutlinedTwoTone'
import Download from '../entry/Download'
import Quantity from '../Quantity'
import FilePreview from './FilePreview'
import { archiveAdaptorFactory } from './ArchiveBrowser'

const FileBrowser = React.memo(({uploadId, path, rootTitle, highlightedItem}) => {
  const adaptor = new RawDirectoryAdaptor(uploadId, path, rootTitle, highlightedItem)
  return <Browser adaptor={adaptor} />
})
FileBrowser.propTypes = {
  uploadId: PropTypes.string.isRequired,
  path: PropTypes.string.isRequired,
  rootTitle: PropTypes.string.isRequired,
  highlightedItem: PropTypes.string
}
export default FileBrowser

class RawDirectoryAdaptor extends Adaptor {
  constructor(uploadId, path, title, highlightedItem) {
    super()
    this.uploadId = uploadId
    this.path = path
    this.title = title
    this.highlightedItem = highlightedItem
    this.data = undefined // Will be set by RawDirectoryContent component when loaded
  }
  async initialize(api) {
    if (this.data === undefined) {
      const encodedPath = this.path.split('/').map(segment => encodeURIComponent(segment)).join('/')
      const response = await api.get(`/uploads/${this.uploadId}/rawdir/${encodedPath}?include_entry_info=true&page_size=500`)
      const elementsByName = {}
      response.directory_metadata.content.forEach(element => { elementsByName[element.name] = element })
      this.data = {response, elementsByName}
    }
  }
  itemAdaptor(key) {
    const ext_path = this.path ? this.path + '/' + key : key
    const element = this.data.elementsByName[key]
    if (element) {
      if (element.is_file) {
        return new RawFileAdaptor(this.uploadId, ext_path, element)
      } else {
        return new RawDirectoryAdaptor(this.uploadId, ext_path, key)
      }
    }
    throw new Error('Bad path: ' + key)
  }
  render() {
    return <RawDirectoryContent uploadId={this.uploadId} path={this.path} title={this.title} highlightedItem={this.highlightedItem}/>
  }
}

function RawDirectoryContent({uploadId, path, title, highlightedItem}) {
  const lane = useContext(laneContext)
  const encodedPath = path.split('/').map(segment => encodeURIComponent(segment)).join('/')

  if (lane.adaptor.data === undefined) {
    return <Content key={path}><Typography>loading ...</Typography></Content>
  } else {
    // Data loaded
    const downloadurl = `uploads/${uploadId}/raw/${encodedPath}?compress=true`
    const segments = path.split('/')
    const lastSegment = segments[segments.length - 1]
    const downloadFilename = `${uploadId}${lastSegment ? ' - ' + lastSegment : ''}.zip`
    return (
      <Content key={path}>
        <Title
          title={title}
          label="folder"
          tooltip={path}
          actions={(
            <Download
              component={IconButton} disabled={false}
              size="small"
              tooltip="download this folder"
              url={downloadurl}
              fileName={downloadFilename}>
              <DownloadIcon />
            </Download>
          )}
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
                <Typography>{element.name}</Typography>
              </Item>
            ))
          }
          {
            lane.adaptor.data.response.pagination.total > 500 &&
              <Typography color="error">Only showing the first 500 rows</Typography>
          }
        </Compartment>
      </Content>)
  }
}
RawDirectoryContent.propTypes = {
  uploadId: PropTypes.string.isRequired,
  path: PropTypes.string.isRequired,
  title: PropTypes.string.isRequired,
  highlightedItem: PropTypes.bool
}

class RawFileAdaptor extends Adaptor {
  constructor(uploadId, path, data) {
    super()
    this.uploadId = uploadId
    this.path = path
    this.data = data
  }
  async itemAdaptor(key, api) {
    if (key === 'archive') {
      if (!this.data.archive) {
        const response = await api.get(`entries/${this.data.entry_id}/archive`)
        this.data.archive = response.data.archive
      }

      return archiveAdaptorFactory(this.data.archive)
    }
  }
  render() {
    return <RawFileContent uploadId={this.uploadId} path={this.path} data={this.data} key={this.path}/>
  }
}

function RawFileContent({uploadId, path, data}) {
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
  const encodedPath = path.split('/').map(segment => encodeURIComponent(segment)).join('/')
  const downloadUrl = `uploads/${uploadId}/raw/${encodedPath}?ignore_mime_type=true`
  return (
    <Content
      key={path}
      display="flex" flexDirection="column" height="100%"
      paddingTop={0} paddingBottom={0} maxWidth={600} minWidth={600}
    >
      <Box paddingTop={1}>
        <Title
          title={data.name}
          label="file"
          tooltip={path}
          actions={(
            <Download component={IconButton} disabled={false}
              size="small"
              tooltip="download this file"
              url={downloadUrl}
              fileName={data.name}>
              <DownloadIcon />
            </Download>
          )}
        />
      </Box>
      <Compartment>
        {/* <Quantity quantity="filename" data={{filename: data.name}} label="file name" noWrap ellipsisFront withClipboard />
        <Quantity quantity="path" data={{path: path}} label="full path" noWrap ellipsisFront withClipboard /> */}
        <Quantity quantity="filesize" data={{filesize: niceSize}} label="file size" />
        { data.parser_name && <Quantity quantity="parser" data={{parser: data.parser_name}} />}
        { data.parser_name && <Quantity quantity="entryId" data={{entryId: data.entry_id}} noWrap withClipboard />}
      </Compartment>
      {data.entry_id && <Compartment>
        <Item itemKey="archive"><Typography>processed data</Typography></Item>
      </Compartment>}
      <Box marginTop={2}/>
      <Box flexGrow={1} overflow="hidden">
        <FilePreview uploadId={uploadId} path={path} size={data.size}/>
      </Box>
      <Box paddingBottom={1}/>
    </Content>)
}
RawFileContent.propTypes = {
  uploadId: PropTypes.string.isRequired,
  path: PropTypes.string.isRequired,
  data: PropTypes.object.isRequired
}
