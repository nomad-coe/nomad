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
import React, { useEffect, useContext, useState } from 'react'
import { Typography, Grid, IconButton } from '@material-ui/core'
import { useErrors } from '../errors'
import Browser, { Item, Content, Adaptor, laneContext } from './Browser'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import Download from '../entry/Download'
import Quantity from '../Quantity'
import { useApi } from '../api'

const FileBrowser = React.memo(({uploadId, path}) => {
  const adaptor = new RawDirectoryAdaptor(uploadId, path)
  return <Browser adaptor={adaptor} />
})
export default FileBrowser

class RawDirectoryAdaptor extends Adaptor {
  constructor(uploadId, path) {
    super()
    this.uploadId = uploadId
    this.path = path
    this.data = undefined  // Will be set by RawDirectoryContent component when loaded
  }
  isLoaded() {
    return this.data !== undefined
  }
  itemAdaptor(key) {
    const ext_path = this.path ? this.path + '/' + key : key
    for(let element of this.data) {
      if(element.name === key) {
        if(element.is_file)
          return new FilePreviewAdaptor(this.uploadId, ext_path, element)
        else
          return new RawDirectoryAdaptor(this.uploadId, ext_path)
      }
    }
    throw new Error('Bad path: ' + key)
  }
  render() {
    return <RawDirectoryContent uploadId={this.uploadId} path={this.path}/>
  }
}

function RawDirectoryContent({uploadId, path}) {
  const {api} = useApi()
  const {raiseError} = useErrors()
  const lane = useContext(laneContext)
  const [loadedPath, setLoadedPath] = useState(null)

  useEffect(() => {
    if (loadedPath !== path) {
      api.get(`/uploads/${uploadId}/raw/` + encodeURIComponent(path))
      .then(response => {
        lane.adaptor.data = response.content
        setLoadedPath(path)
        lane.update()
      })
      .catch(error => {
        raiseError(error)
      })
    }
  }, [uploadId, path, api, lane, loadedPath, setLoadedPath, raiseError])

  if (loadedPath !== path) {
    return <Content key={path}><Typography>loading ...</Typography></Content>
  } else {
    // Data loaded
    return (
      <Content key={path}>
        {lane.adaptor.data.map(element => (
          <Item itemKey={element.name} key={path ? path + '/' + element.name : element.name}>
            {
              element.is_file
              ? <Typography>{element.name}</Typography>
              : <Typography><b>{element.name + '/'}</b></Typography>
            }
          </Item>))}
      </Content>)
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
      // A nicer, human-readable size string
      let size = this.data.size
      let niceSize, unit, factor
      if(size > 1e9) {
        [unit, factor] = ['GB', 1e9]
      } else if (size > 1e6) {
        [unit, factor] = ['MB', 1e6]
      } else if (size > 1e3) {
        [unit, factor] = ['kB', 1e3]
      }
      if(unit) {
        if(size / factor > 100)
          // No decimals
          niceSize = `${Math.round(size / factor)} ${unit} (${size} bytes)`
        else
          // One decimal
          niceSize = `${Math.round(size * 10 / factor) / 10} ${unit} (${size} bytes)`
      } else {
        // Unit = bytes
        niceSize = `${size} bytes`
      }
      let encodedPath = this.path.split('/').map(segment => encodeURIComponent(segment)).join('/')
      let downloadUrl = `uploads/${this.uploadId}/raw/${encodedPath}`

      return (
        <Content key={this.path}>
          <Grid container justifyContent="flex-end">
            <Grid item>
              <Download component={IconButton} disabled={false}
                  color="secondary"
                  tooltip="download this file"
                  url={downloadUrl}
                  fileName={this.data.name}>
                <DownloadIcon />
              </Download>
            </Grid>
          </Grid>
          <Quantity quantity="filename" data={{filename: this.data.name}} label="file name" noWrap ellipsisFront withClipboard />
          <Quantity quantity="path" data={{path: this.path}} label="full path" noWrap ellipsisFront withClipboard />
          <Quantity quantity="filesize" data={{filesize: niceSize}} label="file size" />
        </Content>
      )
    }
}