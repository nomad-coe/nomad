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
import React, { useState, useRef, useCallback, useEffect} from 'react'
import PropTypes from 'prop-types'
import { Typography, makeStyles, Button, Box } from '@material-ui/core'
import { useErrors } from '../errors'
import ReactJson from 'react-json-view'
import InfiniteScroll from 'react-infinite-scroller'
import { useApi } from '../api'
import { apiBase } from '../../config'

const useFilePreviewStyles = makeStyles(theme => ({
  scrollableContainer: {
    width: '100%',
    height: '100%',
    display: 'inline-block',
    overflow: 'auto'
  },
  imgDiv: {
    width: '100%',
    height: '100%',
    position: 'relative'
  },
  imgElement: {
    maxWidth: '100%',
    maxHeight: '100%',
    position: 'absolute',
    top: 0,
    bottom: 0,
    left: 0,
    right: 0,
    margin: 'auto'
  }
}))

/* Viewer definitions */
const viewerText = {
  name: 'text',
  fileExtensions: ['txt', 'yaml', 'yml'],
  maxSizeAutoPreview: 10e6,
  render: ({uploadId, path}) => {
    return <FilePreviewText uploadId={uploadId} path={path}/>
  }
}
const viewerImg = {
  name: 'image',
  fileExtensions: ['png', 'jpg', 'jpeg', 'gif', 'bmp', 'svg'],
  maxSizeAutoPreview: 10e6,
  requiresUrlWithAuth: true,
  render: ({classes, url}) => {
    return <div className={classes.imgDiv}><img src={url} className={classes.imgElement} alt="Loading..."/></div>
  }
}
const viewerJSON = {
  name: 'json',
  fileExtensions: ['json'],
  maxSizeAutoPreview: 10e6,
  requiresLoadedData: true,
  render: ({classes, data}) => {
    if (typeof data.current === 'string') {
      data.current = JSON.parse(data.current)
    }
    return (
      <Box className={classes.scrollableContainer}>
        <ReactJson src={data.current} enableClipboard={false} collapsed={1} displayObjectSize={false}/>
      </Box>
    )
  }
}
const viewers = [viewerText, viewerImg, viewerJSON]

export default function FilePreview({uploadId, path, size}) {
  const classes = useFilePreviewStyles()
  const {api, user} = useApi()
  const {raiseError} = useErrors()

  // Determine viewer to use and if we should preview automatically, based on extension and size
  const fileExtension = path.split('.').pop().toLowerCase()
  let autoPreview = false
  let selectedViewer = viewerText
  for (const viewer of viewers) {
    if (viewer.fileExtensions.includes(fileExtension)) {
      selectedViewer = viewer
      autoPreview = size < viewer.maxSizeAutoPreview
      break
    }
  }

  const [preview, setPreview] = useState(autoPreview)
  const data = useRef()
  const [dataLoaded, setDataLoaded] = useState(false)

  const encodedPath = path.split('/').map(segment => encodeURIComponent(segment)).join('/')
  let fullUrl = `${apiBase}/v1/uploads/${uploadId}/raw/${encodedPath}`
  if (fullUrl.startsWith('/')) {
    fullUrl = `${window.location.origin}${fullUrl}`
  }
  const [fullUrlWithAuth, setFullUrlWithAuth] = useState(undefined)

  useEffect(() => {
    if (preview && user && !fullUrlWithAuth && selectedViewer.requiresUrlWithAuth) {
      // Need to fetch signature token for the viewer
      api.get('/auth/signature_token')
        .then(response => {
          const fullUrlWithAuth = new URL(fullUrl)
          fullUrlWithAuth.searchParams.append('signature_token', response.signature_token)
          setFullUrlWithAuth(fullUrlWithAuth.href)
        })
        .catch(raiseError)
    }
    if (preview && selectedViewer.requiresLoadedData && data.current === undefined) {
      // Need to load the file data content for the viewer
      data.current = null
      api.get(fullUrl)
        .then(response => {
          data.current = response
          setDataLoaded(true)
        })
        .catch(raiseError)
    }
  }, [preview, selectedViewer, user, fullUrl, fullUrlWithAuth, setFullUrlWithAuth, data, dataLoaded, setDataLoaded, api, raiseError])

  if (!preview) {
    return (
      <Box margin={2} textAlign="center">
        <Button onClick={() => setPreview(true)} variant="contained" size="small" color="primary">
          Preview with {selectedViewer.name} viewer
        </Button>
      </Box>
    )
  }

  const url = user ? fullUrlWithAuth : fullUrl
  if ((selectedViewer.requiresUrlWithAuth && !url) || (selectedViewer.requiresLoadedData && !dataLoaded)) {
    // Not ready to invoke the viewer yet
    return <Typography>Loading...</Typography>
  }

  try {
    return selectedViewer.render({uploadId, path, url, data, classes})
  } catch (error) {
    return <Typography>Failed to open viewer</Typography>
  }
}
FilePreview.propTypes = {
  uploadId: PropTypes.string.isRequired,
  path: PropTypes.string.isRequired,
  size: PropTypes.number.isRequired
}

const useFilePreviewTextStyles = makeStyles(theme => ({
  containerDiv: {
    width: '100%',
    height: '100%',
    overflow: 'auto',
    backgroundColor: theme.palette.primary.dark
  },
  fileContents: {
    margin: 0,
    padding: 0,
    display: 'inline-block',
    color: theme.palette.primary.contrastText,
    fontFamily: 'Consolas, "Liberation Mono", Menlo, Courier, monospace',
    fontSize: 12,
    minWidth: '100%'
  }
}))
function FilePreviewText({uploadId, path}) {
  const classes = useFilePreviewTextStyles()
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [contents, setContents] = useState(null)
  const [hasMore, setHasMore] = useState(true)
  const [loading, setLoading] = useState(false)
  const containerRef = React.createRef()

  const encodedPath = path.split('/').map(segment => encodeURIComponent(segment)).join('/')

  const loadMore = useCallback(() => {
    // The infinite scroll component has the issue if calling load more whenever it
    // gets updates, therefore calling this infinitely before it gets any chances of
    // receiving the results (https://github.com/CassetteRocks/react-infinite-scroller/issues/163).
    // Therefore, we have to set hasMore to false first and set it to true again after
    // receiving actual results.    setLoading(true)
    if (hasMore && !loading) {
      api.get(
        `/uploads/${uploadId}/raw/${encodedPath}`,
        {
          offset: contents?.length || 0,
          length: 16 * 1024,
          decompress: true,
          ignore_mime_type: true
        },
        {transformResponse: []})
        .then(contents => {
          setContents(old => (old || '') + (contents || ''))
          setHasMore(contents?.length === 16 * 1024)
        })
        .catch(raiseError)
        .finally(() => setLoading(false))
    }
  }, [uploadId, encodedPath, loading, hasMore, setHasMore, setContents, api, raiseError, contents])

  if (loading && !contents) {
    return <Typography>Loading ...</Typography>
  }
  return (
    <div ref={containerRef} className={classes.containerDiv}>
      <InfiniteScroll
        pageStart={0}
        loadMore={loadMore}
        hasMore={hasMore}
        useWindow={false}
        getScrollParent={() => containerRef.current}
      >
        <pre className={classes.fileContents}>
          {contents}
          &nbsp;
        </pre>
      </InfiniteScroll>
    </div>
  )
}
FilePreviewText.propTypes = {
  uploadId: PropTypes.string.isRequired,
  path: PropTypes.string.isRequired
}
