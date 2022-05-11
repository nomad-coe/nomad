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
import { Box, Button, Dialog, DialogActions, DialogContent, DialogTitle, Typography } from '@material-ui/core'
import React, { useState } from 'react'
import PropTypes from 'prop-types'
import { useApi } from '../api'
import { useErrors } from '../errors'
import exampleUploads from '../../exampleUploads.json'
import Markdown from '../Markdown'
import { useHistory } from 'react-router-dom'

const ExampleUploadDialog = React.memo(function ExampleUploadDialog(props) {
  const { onSelect, ...dialogProps } = props
  return (
    <Dialog {...dialogProps}>
      <DialogTitle id="Uploads-dialogBox">Select a sample Upload</DialogTitle>
      <DialogContent>
        {Object.keys(exampleUploads).map(key => {
          const exampleUpload = exampleUploads[key]
          return (
            <Box key={key} display="flex" flexDirection="row" marginBottom={2} alignItems="baseline">
              <Box marginRight={2} flexGrow={1}>
                <Typography variant="h6">{exampleUpload.title}</Typography>
                <Markdown>{`${exampleUpload.description}`}</Markdown>
              </Box>
              <Button
                variant="contained" color="primary" size="small"
                onClick={() => onSelect(key)}
              >
                add
              </Button>
            </Box>
          )
        })}
      </DialogContent>
      <DialogActions>
        <Button onClick={dialogProps.onClose}>cancel</Button>
      </DialogActions>
    </Dialog>
  )
})

ExampleUploadDialog.propTypes = {
  onClose: PropTypes.func.isRequired,
  onSelect: PropTypes.func.isRequired
}

export default function ExampleUploadButton(props) {
  const {api} = useApi()
  const {raiseError} = useErrors()
  const history = useHistory()

  const [isOpen, setOpen] = useState(false)

  const handleClickOpen = () => {
    setOpen(true)
  }

  const handleClose = () => setOpen(false)

  const handleSelect = (value) => {
    const exampleUpload = exampleUploads[value]
    api.post(`/uploads?local_path=${exampleUpload.path}&upload_name=${exampleUpload.title}`)
      .then((data) => {
        history.push(`/user/uploads/upload/id/${data.upload_id}`)
      })
      .catch(raiseError)
      .finally(() => {
        setOpen(false)
      })
  }

  return (
    <div>
      <Button variant="contained" onClick={handleClickOpen}>
        Add Example Uploads
      </Button>
      <ExampleUploadDialog
        open={isOpen}
        onClose={handleClose}
        onSelect={handleSelect}
      />
    </div>
  )
}
ExampleUploadButton.propTypes = {
  isDisable: PropTypes.bool
}
