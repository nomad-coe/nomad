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

import React, { useRef, useState } from 'react'
import PropTypes from 'prop-types'
import { Button, Dialog, DialogActions, DialogContent, DialogTitle, makeStyles } from '@material-ui/core'
import InfiniteScroll from 'react-infinite-scroller'

const useFilePreviewStyles = makeStyles(theme => ({
  contents: {
    width: '85%',
    overflowX: 'auto',
    color: theme.palette.primary.contrastText,
    backgroundColor: theme.palette.primary.dark,
    marginTop: theme.spacing(1),
    padding: '3px 6px',
    fontFamily: 'Consolas, "Liberation Mono", Menlo, Courier, monospace',
    fontSize: 12
  }
}))

export default function FilePreview({uploadId, path, onClose}) {
  const scrollRef = useRef()
  const [content] = useState(null)
  const classes = useFilePreviewStyles()

  return <Dialog open>
    <DialogTitle>File preview: {path}</DialogTitle>
    <DialogContent ref={scrollRef}>
      <InfiniteScroll
        className={classes.fileContents}
        pageStart={0}
        loadMore={() => console.log('### load more')}
        hasMore={false}
        useWindow={false}
        getScrollParent={() => scrollRef.current}
      >
        <pre style={{margin: 0}}>
          {`${content}`}
          &nbsp;
        </pre>
      </InfiniteScroll>
    </DialogContent>
    <DialogActions>
      <Button onClick={onClose} color="primary">
        Cancel
      </Button>
    </DialogActions>
  </Dialog>
}

FilePreview.propTypes = {
  uploadId: PropTypes.string.isRequired,
  path: PropTypes.string.isRequired,
  onClose: PropTypes.func
}
