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
import React, {useState, useCallback, useEffect} from 'react'
import PropTypes from 'prop-types'
import {
  makeStyles,
  useTheme,
  FormGroup,
  FormControlLabel,
  Checkbox,
  FormLabel,
  IconButton,
  Divider,
  Typography,
  Tooltip
} from '@material-ui/core'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import Download from './Download'
import ReloadIcon from '@material-ui/icons/Cached'
import ViewIcon from '@material-ui/icons/Search'
import InfiniteScroll from 'react-infinite-scroller'
import { ScrollContext } from '../nav/Navigation'
import { useApi } from '../api'
import { useErrors } from '../errors'

const useStyles = makeStyles(theme => ({
  root: {},
  formLabel: {
    padding: theme.spacing(2)
  },
  shownFile: {
    color: theme.palette.secondary.main,
    overflowX: 'hidden',
    textOverflow: 'ellipsis',
    whiteSpace: 'nowrap',
    direction: 'rtl',
    textAlign: 'left'
  },
  fileContents: {
    width: '85%',
    overflowX: 'auto',
    color: theme.palette.primary.contrastText,
    backgroundColor: theme.palette.primary.dark,
    marginTop: theme.spacing(1),
    padding: '3px 6px',
    fontFamily: 'Consolas, "Liberation Mono", Menlo, Courier, monospace',
    fontSize: 12
  },
  fileError: {
    marginTop: 16,
    padding: 8
  },
  fileNameFormGroup: {
    display: 'flex',
    flexWrap: 'nowrap'
  },
  fileNameFormGroupLabel: {
    flexGrow: 1,
    overflowX: 'hidden',
    marginRight: 0
  },
  fileNameLabel: {
    overflowX: 'hidden',
    textOverflow: 'ellipsis',
    whiteSpace: 'nowrap',
    direction: 'rtl',
    textAlign: 'left'
  },
  mainfileLabel: {
    fontSize: '10pt',
    fontWeight: 'bold'
  }
}))

function label(file) {
  return file.split('/').reverse()[0]
}

export default function RawFiles({data, entryId}) {
  const theme = useTheme()
  const classes = useStyles(theme)
  const [selectedFiles, setSelectedFiles] = useState([])
  const [shownFile, setShownFile] = useState(null)
  const [fileContents, setFileContents] = useState(null)
  const [files, setFiles] = useState(null)
  const [loading, setLoading] = useState(false)
  const [doesNotExist, setDoesNotExist] = useState(false)
  const {raiseError} = useErrors()
  const {api, user} = useApi()

  useEffect(() => {
    setSelectedFiles([])
    setShownFile(null)
    setFileContents(null)
    setFiles(null)
    setLoading(false)
    setDoesNotExist(false)
  }, [api, entryId])

  const update = useCallback(() => {
    // this might accidentally happen, when the user logs out and the ids aren't
    // necessarily available anymore, but the component is still mounted
    if (!entryId) {
      return
    }

    api.getRawFileListFromEntry(entryId).then(data => {
      const files = data.data.files.map(file => file.path)
      if (files.length > 500) {
        raiseError('There are more than 500 files in this entry. We can only show the first 500.')
      }
      setFiles(files)
      setLoading(false)
    }).catch(error => {
      setFiles(null)
      setLoading(false)
      if (error.name === 'DoesNotExist') {
        setDoesNotExist(true)
      } else {
        raiseError(error)
      }
    })
    setLoading(true)
  }, [entryId, raiseError, api])

  const handleSelectFile = useCallback((file) => {
    setSelectedFiles(prevState => {
      const index = prevState.indexOf(file)
      if (index === -1) {
        return [file, ...prevState]
      } else {
        return [...prevState.slice(0, index), ...prevState.slice(index + 1)]
      }
    })
  }, [])

  const handleFileClicked = useCallback(file => {
    setShownFile(file)
    setFileContents(null)
    api.get(
      `/entries/${entryId}/raw/${file.split('/').reverse()[0]}`,
      {length: 16 * 1024, decompress: true},
      {transformResponse: []})
      .then(contents => setFileContents({
        hasMore: true,
        contents: contents
      }))
      .catch(raiseError)
  }, [api, raiseError, entryId])

  const handleLoadMore = useCallback((page) => {
    // The infinite scroll component has the issue if calling load more whenever it
    // gets updates, therefore calling this infinitely before it gets any chances of
    // receiving the results (https://github.com/CassetteRocks/react-infinite-scroller/issues/163).
    // Therefore, we have to set hasMore to false first and set it to true again after
    // receiving actual results.
    setFileContents(prevState => { return {...prevState, hasMore: false} })
    const initialEntryId = entryId

    if (fileContents.contents.length < (page + 1) * 16 * 1024) {
      api.get(
        `/entries/${entryId}/raw/${shownFile.split('/').reverse()[0]}`,
        {offset: page * 16 * 1024, length: 16 * 1024, decompress: true},
        {transformResponse: []})
        .then(contents => {
          // The back-button navigation might cause a scroll event, might cause to loadmore,
          // will set this state, after navigation back to this page, but potentially
          // different entry.
          if (initialEntryId === entryId) {
            setFileContents({
              hasMore: contents.length > 0,
              contents: ((fileContents && fileContents.contents) || '') + contents
            })
          }
        })
        .catch(error => {
          setFileContents(null)
          setShownFile(null)
          raiseError(error)
        })
    }
  }, [api, entryId, shownFile, fileContents, raiseError])

  const filterPotcar = useCallback((file) => {
    if (file.substring(file.lastIndexOf('/')).includes('POTCAR') && !file.endsWith('.stripped')) {
      return user && data.main_author.user_id === user.sub
    } else {
      return true
    }
  }, [data, user])

  const availableFiles = files || data.files || []
  const someSelected = selectedFiles.length > 0
  const allSelected = availableFiles.length === selectedFiles.length && someSelected

  if (doesNotExist) {
    return <Typography>
      The uploaded raw files for this entry do not exist. This is most likely a NOMAD
      issue. Please inform us, if this error persists.
    </Typography>
  }

  const file = path => path.substring(path.lastIndexOf('/') + 1)

  let downloadUrl
  if (selectedFiles.length === 1) {
    // download the individual file
    downloadUrl = `entries/${entryId}/raw/${file(selectedFiles[0])}`
  } else if (selectedFiles.length === availableFiles.length) {
    // use an endpoint that downloads all files of the entry
    downloadUrl = `entries/${entryId}/raw`
  } else if (selectedFiles.length > 0) {
    // download specific files
    const query = selectedFiles.map(file).map(f => `include_files=${encodeURIComponent(f)}`).join('&')
    downloadUrl = `entries/${entryId}/raw?${query}`
  }

  return (
    <div className={classes.root}>
      <FormGroup row>
        <FormControlLabel
          label="select all" style={{flexGrow: 1}}
          control={
            <Checkbox value="select_all" checked={allSelected}
              indeterminate={!allSelected && someSelected}
              onChange={() => setSelectedFiles(allSelected ? [] : availableFiles.slice())}
            />
          }
        />
        {!files
          ? <Tooltip title="check for more files">
            <IconButton onClick={() => update()}>
              <ReloadIcon />
            </IconButton>
          </Tooltip> : ''
        }
        <FormLabel className={classes.formLabel}>
          {selectedFiles.length}/{availableFiles.length} files selected
        </FormLabel>
        <Download component={IconButton} disabled={selectedFiles.length === 0}
          color="secondary"
          tooltip="download selected files"
          url={downloadUrl}
          fileName={selectedFiles.length === 1 ? label(selectedFiles[0]) : `${entryId}.zip`}
        >
          <DownloadIcon />
        </Download>
      </FormGroup>
      <Divider />
      <div style={{display: 'flex', flexDirection: 'row'}}>
        <div style={{width: '25%'}}>
          {availableFiles.filter(filterPotcar).map((file, index) => (
            <FormGroup row key={index} className={classes.fileNameFormGroup}>
              <Tooltip title={file + (index === 0 ? ', the main (output) file of this entry' : '')}>
                <FormControlLabel
                  style={{flexGrow: 1, overflowX: 'hidden', textOverflow: 'ellipsis'}}
                  label={<span>{label(file)}{index === 0 ? <b className={classes.mainfileLabel}> (mainfile)</b> : ''}</span>}
                  classes={{
                    root: classes.fileNameFormGroupLabel,
                    label: file === shownFile ? classes.shownFile : classes.fileNameLabel}}
                  control={
                    <Checkbox
                      disabled={loading}
                      checked={selectedFiles.indexOf(file) !== -1}
                      onChange={() => handleSelectFile(file)} value={file}
                    />}
                />
              </Tooltip>
              <Tooltip title='Show contents'>
                <IconButton onClick={() => handleFileClicked(file)} color={file === shownFile ? 'secondary' : 'default'}>
                  <ViewIcon />
                </IconButton>
              </Tooltip>
            </FormGroup>
          ))}
        </div>
        {fileContents && fileContents.contents !== null &&
          <ScrollContext.Consumer>
            {scroll =>
              <InfiniteScroll
                className={classes.fileContents}
                pageStart={0}
                loadMore={handleLoadMore}
                hasMore={fileContents.hasMore}
                useWindow={false}
                getScrollParent={() => scroll.scrollParentRef.current}
              >
                <pre style={{margin: 0}}>
                  {`${fileContents.contents}`}
                  &nbsp;
                </pre>
              </InfiniteScroll>
            }
          </ScrollContext.Consumer>
        }
        {fileContents && fileContents.contents === null &&
          <div className={classes.fileError}>
            <Typography color="error">
              Cannot display file due to unsupported file format.
            </Typography>
          </div>
        }
      </div>
    </div>
  )
}
RawFiles.propTypes = {
  entryId: PropTypes.string.isRequired,
  data: PropTypes.object
}
