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
import React, { useState, useRef, useCallback, useEffect, useLayoutEffect} from 'react'
import PropTypes from 'prop-types'
import { Typography, makeStyles, Button, Box } from '@material-ui/core'
import { useErrors } from '../errors'
import ReactJson from 'react-json-view'
import { Document, Page, pdfjs } from 'react-pdf'
import InfiniteScroll from 'react-infinite-scroller'
import { useApi } from '../api'
import { apiBase } from '../../config'
import { parseCifStructures } from 'crystcif-parse'
import H5Web from '../visualization/H5Web'
import Structure from '../visualization/Structure'
import { isWaitingForUpdateTestId } from '../../utils'
import {useLane} from './Browser'

pdfjs.GlobalWorkerOptions.workerSrc = `https://cdnjs.cloudflare.com/ajax/libs/pdf.js/${pdfjs.version}/pdf.worker.js`

const useFilePreviewStyles = makeStyles(theme => ({
  scrollableContainer: {
    width: '100%',
    height: '100%',
    display: 'inline-block',
    overflow: 'auto'
  },
  img: {
    maxHeight: '100%',
    maxWidth: '100%'
  }
}))

/* Viewer definitions */
const viewerText = {
  name: 'text',
  fileExtensions: ['txt', 'yaml', 'yml', 'csv', 'xml', 'dat'],
  maxSizePreview: 1e10, // Effectively infinite
  maxSizeAutoPreview: 1e10, // Effectively infinite
  width: 700,
  render: ({uploadId, path}) => {
    return <FilePreviewText uploadId={uploadId} path={path}/>
  }
}
const viewerImg = {
  name: 'image',
  fileExtensions: ['png', 'jpg', 'jpeg', 'gif', 'bmp', 'svg'],
  maxSizeAutoPreview: 10e6,
  render: ({classes, url, setFailedToPreview}) => {
    return <img
      src={url} className={classes.img}
      alt="Loading..." onError={() => setFailedToPreview(true)}
    />
  }
}
const viewerJSON = {
  name: 'json',
  fileExtensions: ['json'],
  maxSizeAutoPreview: 10e6,
  requiresLoadedData: true,
  width: 'fit-content',
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
const viewerPDF = {
  name: 'pdf',
  fileExtensions: ['pdf'],
  maxSizeAutoPreview: 10e6,
  width: 700,
  render: ({url, setFailedToPreview}) => {
    return <FilePreviewPdf
      file={{url: url, withCredentials: true}}
      error={(error) => {
        console.log(error)
        setFailedToPreview(true)
      }}
    />
  }
}

const viewerHdfFive = {
  name: 'hdf5',
  fileExtensions: ['h4', 'hd4', 'hdf4', 'hdf', 'h5', 'hd5', 'hdf5', 'hd5', 'nxs', 'h5oina', 'edaxh5', 'emd', 'dream3d'],
  maxSizeAutoPreview: 10e6,
  width: 'fit-content',
  render: ({uploadId, path}) => <H5Web upload_id={uploadId} filename={path}/>
}
const viewerCif = {
  name: 'cif',
  fileExtensions: ['cif'],
  maxSizeAutoPreview: 1e5,
  requiresLoadedData: true,
  width: 500,
  render: ({data}) => {
    if (typeof data.current === 'string') {
      const cifData = parseCifStructures(data.current)
      const cifStructure = cifData[Object.keys(cifData)[0]]
      const structureData = {
        cell: cifStructure.get_cell(),
        pbc: cifStructure.get_pbc(),
        positions: cifStructure.get_positions(),
        species: cifStructure.get_chemical_symbols()
      }
      data.current = structureData
    }
    return (
      <div style={{height: 500, width: 500}}>
        <Structure data={data.current} />
      </div>
    )
  }
}

export const viewers = [viewerText, viewerImg, viewerJSON, viewerPDF, viewerHdfFive, viewerCif]

const FilePreview = React.memo(({uploadId, path, size}) => {
  const classes = useFilePreviewStyles()
  const {api, user} = useApi()
  const {raiseError} = useErrors()

  // Determine viewer to use and if we should preview automatically, based on extension and size
  const fileExtension = path.split('.').pop().toLowerCase()
  let autoPreview = true
  let selectedViewer = viewerText
  for (const viewer of viewers) {
    if (viewer.fileExtensions.includes(fileExtension)) {
      if (size < (selectedViewer.maxSizePreview || 50e6)) {
        selectedViewer = viewer
        autoPreview = size < viewer.maxSizeAutoPreview
        break
      }
    }
  }

  const [preview, setPreview] = useState(autoPreview)
  const [failedToPreview, setFailedToPreview] = useState(false)
  const [useFallbackViewer, setUseFallbackViewer] = useState(false)

  const data = useRef()
  const [dataLoaded, setDataLoaded] = useState(false)

  const encodedPath = path.split('/').map(segment => encodeURIComponent(segment)).join('/')
  let url = `${apiBase}/v1/uploads/${uploadId}/raw/${encodedPath}`
  if (url.startsWith('/')) {
    url = `${window.location.origin}${url}`
  }

  useEffect(() => {
    if (preview && !failedToPreview && selectedViewer.requiresLoadedData && data.current === undefined) {
      // Need to load the file data content for the viewer
      data.current = null
      api.get(url)
        .then(response => {
          data.current = response
          setDataLoaded(true)
        })
        .catch(error => {
          setFailedToPreview(true)
          raiseError(error)
        })
    }
  }, [preview, failedToPreview, selectedViewer, user, url, data, dataLoaded, setDataLoaded, api, raiseError])

  if (!preview) {
    return (
      <Box margin={2} textAlign="center">
        <Button onClick={() => setPreview(true)} variant="contained" size="small" color="primary">
          Preview with {selectedViewer.name} viewer
        </Button>
      </Box>
    )
  }

  let content
  if (selectedViewer.requiresLoadedData && !dataLoaded) {
    // Not ready to invoke the viewer yet
    content = <Typography>Loading...</Typography>
  } else if (!failedToPreview) {
    try {
      content = selectedViewer.render({uploadId, path, url, data, setFailedToPreview, classes})
    } catch (error) {
      // TODO
    }
  } else if (!useFallbackViewer) {
    // Selected viewer failed
    content = (
      <Box textAlign="center">
        <Typography color="error">Failed to open with {selectedViewer.name} viewer. Bad file format?</Typography>
        <Button onClick={() => setUseFallbackViewer(true)} variant="contained" size="small" color="primary">
          Open with text viewer
        </Button>
      </Box>
    )
  } else {
    // Use the text viewer as last resort
    content = viewerText.render({uploadId, path, url, data, setFailedToPreview, classes})
  }
  return <Box marginBottom={1} flexGrow={1} width={`${selectedViewer.width || 500}px`} overflow="hidden">
    {content}
  </Box>
})
FilePreview.propTypes = {
  uploadId: PropTypes.string.isRequired,
  path: PropTypes.string.isRequired,
  size: PropTypes.number.isRequired
}
export default FilePreview

const useFilePreviewTextStyles = makeStyles(theme => ({
  containerDiv: {
    height: '100%',
    backgroundColor: theme.palette.primary.dark
  },
  fileContents: {
    margin: theme.spacing(1),
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
  const loading = useRef(false)
  const [loadFailed, setLoadFailed] = useState(false)
  const lane = useLane()

  const encodedPath = path.split('/').map(segment => encodeURIComponent(segment)).join('/')

  const loadMore = useCallback(async () => {
    // The infinite scroll component normally calls the loadMore hook whenever it
    // gets updates, therefore calling this infinitely before it gets any chances of
    // receiving the results (https://github.com/CassetteRocks/react-infinite-scroller/issues/163).
    // Therefore, we use `loading` to keep track of if we are already running a load request.
    if (hasMore && !loading.current && !loadFailed) {
      loading.current = true
      api.get(
        `/uploads/${uploadId}/raw/${encodedPath}`,
        {
          offset: contents?.length || 0,
          length: 16 * 1024,
          decompress: true,
          ignore_mime_type: true
        },
        {
          transformResponse: [],
          responseType: 'blob'
        })
        .then(content => {
          // Note that the length of the string is not necessarily the number of
          // bytes that were transferred. This is why we need to request a blob,
          // get it's size, and then in a second step decode it as string for
          // displaying.
          const nBytes = content.size
          setHasMore(nBytes === 16 * 1024)
          return content.text()
        })
        .then(text => {
          setContents(old => (old || '') + (text || ''))
        })
        .catch(error => {
          setLoadFailed(true)
          setContents(old => (old || '') + '\n\nERROR - COULD NOT READ FILE!')
          raiseError(error)
        })
        .finally(() => { loading.current = false })
    }
  }, [uploadId, encodedPath, loading, loadFailed, setLoadFailed, hasMore, setHasMore, setContents, api, raiseError, contents])

  useEffect(() => {
    // Trigger loading the first chunk
    if (contents === null && !loading.current) {
      loadMore()
    }
  }, [contents, loading, loadMore])

  return (
    <div className={classes.containerDiv}>
      <InfiniteScroll
        pageStart={0}
        loadMore={loadMore}
        hasMore={hasMore}
        useWindow={false}
        getScrollParent={() => lane.containerRef.current}
        {...(contents === null ? {'data-testid': isWaitingForUpdateTestId} : {})}
      >
        <pre className={classes.fileContents}>
          {contents || ''}
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

const useFilePreviewPdfStyles = makeStyles(theme => ({
  containerDiv: {
    height: '100%',
    overflowX: 'hidden',
    overflowY: 'scroll',
    border: '1px solid rgba(0, 0, 0, 0.42)',
    boxSizing: 'border-box'
  },
  pageDiv: {
    border: '5px solid gray'
  }
}))
const FilePreviewPdf = React.memo(props => {
  const classes = useFilePreviewPdfStyles()
  const containerRef = useRef()
  const [numPages, setNumPages] = useState(null)
  const [pageWidth, setPageWidth] = useState()

  useLayoutEffect(() => {
    if (containerRef.current) {
      setPageWidth(containerRef.current.clientWidth - 10) // The -10 is because of the borders
    }
  }, [])

  function onDocumentLoadSuccess({ numPages }) {
    setNumPages(numPages)
  }

  return (
    <div className={classes.containerDiv} ref={containerRef}>
      <Document onLoadSuccess={onDocumentLoadSuccess} renderMode="svg" {...props}>
        {numPages && pageWidth &&
          Array(numPages).fill().map((_, i) =>
            <div className={classes.pageDiv} key={`pdfPageDiv${i + 1}`}>
              <Page pageNumber={i + 1} key={`pdfPage${i + 1}`} width={pageWidth}/>
            </div>)
        }
      </Document>
    </div>)
})
