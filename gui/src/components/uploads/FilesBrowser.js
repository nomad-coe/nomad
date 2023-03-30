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

import React, { useEffect, useMemo, useRef, useState } from 'react'
import PropTypes from 'prop-types'
import { Box, Chip, CircularProgress, Collapse, IconButton, makeStyles, Paper, Tooltip, Typography } from '@material-ui/core'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ChevronRightIcon from '@material-ui/icons/ChevronRight'
import NavigateIcon from '@material-ui/icons/ArrowForward'
import { useApi } from '../api'
import { useErrors } from '../errors'
import FilePreview from './FilePreview'
import { EntryButton } from '../nav/Routes'

const useFolderStyles = makeStyles(theme => ({
  root: {},
  item: {
    '&:hover': {
      backgroundColor: theme.palette.grey[100]
    },
    width: '100%',
    display: 'flex',
    flexDirection: 'raw',
    alignItems: 'center',
    flexWrap: 'nowrap',
    height: 32,
    paddingRight: theme.spacing(1)
  },
  icon: {
  },
  name: {

  },
  tags: {
    // textTransform: 'uppercase',
    marginLeft: theme.spacing(1)
  },
  info: {
    marginLeft: theme.spacing(1)
  },
  actions: {
    marginLeft: theme.spacing(1),
    display: 'flex',
    flexDirection: 'row',
    justifyContent: 'flex-end',
    flexWrap: 'nowrap'
  },
  children: {
    marginLeft: theme.spacing(3)
  }
}))

function FileOrFolder({onToggle, open, hasChildren, children, name, parser, info, entry_id, path}) {
  const classes = useFolderStyles()
  const handleToggle = event => {
    event.stopPropagation()
    if (onToggle) {
      onToggle()
    }
  }
  const iconProps = {
    className: classes.icon,
    fontSize: 'small'
  }
  const icon = hasChildren
    ? (open ? <ExpandMoreIcon {...iconProps} /> : <ChevronRightIcon {...iconProps}/>)
    : <div className={classes.icon} />

  return <div className={classes.root}>
    <div onClick={handleToggle} className={classes.item}>
      {icon}
      <Box marginLeft={icon ? 1 : 0}>
        <Typography className={classes.name}>{name || '/'}</Typography>
      </Box>
      {(info || parser) && <div className={classes.info}>
        {info}
      </div>}
      <div className={classes.tags}>
        {parser && <Chip size="small" label={parser.replace('archive', 'nomad')} color="default" />}
      </div>
      {entry_id && (
        <React.Fragment>
          <Box flexGrow={1} />
          <Tooltip title="Go to the entry page">
            <EntryButton entryId={entry_id} component={IconButton} size="small">
              <NavigateIcon />
            </EntryButton>
          </Tooltip>
        </React.Fragment>
      )}
    </div>
    <Collapse in={open} className={classes.children}>
      {children || 'loading ...'}
    </Collapse>
  </div>
}

FileOrFolder.propTypes = {
  onToggle: PropTypes.func,
  open: PropTypes.bool,
  hasChildren: PropTypes.bool,
  children: PropTypes.arrayOf(PropTypes.object),
  name: PropTypes.string,
  path: PropTypes.string,
  entry_id: PropTypes.string,
  parser: PropTypes.string,
  info: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
}

const useStyles = makeStyles(theme => ({
  root: {
    padding: theme.spacing(2),
    position: 'relative'
  },
  disabled: {
    position: 'absolute',
    left: 0,
    right: 0,
    top: 0,
    bottom: 0,
    borderRadius: 4,
    background: 'rgba(200, 200, 200, 0.5)',
    zIndex: 10
  },
  disabledProgress: {
    position: 'absolute',
    margin: 0,
    top: '50%',
    left: '50%',
    transform: 'translate(-50%, -50%)'
  }
}))

export default function FilesBrower({uploadId, disabled}) {
  const classes = useStyles()
  const {api} = useApi()
  const errors = useErrors()
  const [renderCounter, setRenderCounter] = useState(-1)
  const [previewPath, setPreviewPath] = useState(null)
  const allData = useRef({})

  const fetchData = useMemo(() => (path, open) => {
    async function fetchData() {
      const encodedPath = path.split('/').map(segment => encodeURIComponent(segment)).join('/')
      const results = await api.get(`/uploads/${uploadId}/rawdir/${encodedPath}?page_size=500&include_entry_info=true`)
      allData.current[path] = {
        open: open,
        ...results
      }
      const resultsByPath = {}
      results.directory_metadata.content
        .filter(item => item.is_file)
        .forEach(item => {
          resultsByPath[`${path}/${item.name}`] = item
          if (item.parser_name) {
            item.parser = item.parser_name.replace('parsers/', '').replace('parser/', '')
          }
        })
      setRenderCounter(renderCounter => renderCounter + 1)
    }
    if (!disabled) {
      fetchData().catch(errors.raiseError)
    }
  }, [uploadId, api, errors, setRenderCounter, disabled])

  useEffect(() => {
    if (!disabled) {
      fetchData('')
    }
  }, [fetchData, disabled])

  const handleToggle = (path) => {
    const folderData = allData.current[path]
    if (folderData) {
      // TODO this can be avoided, if this component gets notified about possible file
      // changes behind
      // the api.
      if (!folderData.open) {
        fetchData(path, true)
      }
      folderData.open = !folderData.open
      setRenderCounter(renderCounter + 1)
    } else {
      fetchData(path, true)
    }
  }

  function renderFileOrFolder(path, item) {
    const {is_file} = item
    const data = allData.current[path]
    const pathPrefix = path ? `${path}/` : ''
    const mapContent = item => renderFileOrFolder(`${pathPrefix}${item.name}`, item)
    const props = {
      key: path,
      path: path,
      uploadId: uploadId,
      hasChildren: !is_file,
      open: data?.open,
      children: data?.directory_metadata?.content?.map(mapContent),
      onToggle: is_file ? null : () => handleToggle(path),
      // TODO
      // info: !is_file && data?.content?.length === 0 && <Typography variant="caption">
      //   {'empty directory'}
      // </Typography>,
      ...item
    }

    return <FileOrFolder {...props} />
  }

  const root = allData.current['']

  return <Paper className={classes.root} elevation={disabled ? 0 : undefined}>
    {disabled && <div className={classes.disabled} />}
    {disabled && <div className={classes.disabledProgress}><CircularProgress color="inherit" /></div>}
    {root && renderFileOrFolder('', root)}
    {previewPath && <FilePreview uploadId={uploadId} path={previewPath} onClose={() => setPreviewPath(null)} />}
  </Paper>
}

FilesBrower.propTypes = {
  uploadId: PropTypes.string.isRequired,
  disabled: PropTypes.bool
}
