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
import { Button, Dialog, DialogTitle, makeStyles } from '@material-ui/core'
import React, { useState } from 'react'
import PropTypes from 'prop-types'
import List from '@material-ui/core/List'
import ListItem from '@material-ui/core/ListItem'
import CloudUploadIcon from '@material-ui/icons/CloudUpload'
import { useApi } from '../api'
import { useErrors } from '../errors'
import exampleUploads from '../../ExampleUploads.json'
import Markdown from '../Markdown'

const useStyles = makeStyles(theme => ({
  markDownHeader: {
    color: '#000000',
    fontSize: '18px'
  },
  markDownDescription: {
    color: '#808080',
    fontSize: '15px'
  },
  listItem: {
    flexDirection: 'column',
    maxWidth: '500px',
    alignItems: 'start'
  },
  mainDiv: {
    display: 'flex',
    alignItems: 'center'
  },
  selectButton: {
    marginRight: '5px',
    height: '30px'
  }
}))

const ExampleUploadDialog = React.memo(function ExampleUploadDialog(props) {
  const { exampleUploadData, onClose, openDialog } = props
  const classes = useStyles()

  return (
    <Dialog
      onClose={() => onClose()}
      aria-labelledby="Uploads-dialogBox"
      open={openDialog}
    >
      <DialogTitle id="Uploads-dialogBox">Select a sample Upload</DialogTitle>
      <List>
        {Object.keys(exampleUploadData).map((exampleType) => exampleUploadData[exampleType].map((uploadEntry, i) =>
          <div key={i} className={classes.mainDiv}>
            <ListItem className={classes.listItem}>
              <Markdown className={classes.markDownHeader}>
                {`
                  ${uploadEntry.title}  (${exampleType})
                `}
              </Markdown>
              <Markdown className={classes.markDownDescription}>
                {`
                  ${uploadEntry.description}
                `}
              </Markdown>
            </ListItem>
            <Button
              variant="contained"
              color="primary"
              startIcon={<CloudUploadIcon />}
              onClick={() => onClose(uploadEntry.filePath)}
              className={classes.selectButton}
            >
              Select
            </Button>
          </div>
        ))}
      </List>
    </Dialog>
  )
})

ExampleUploadDialog.propTypes = {
  onClose: PropTypes.func.isRequired,
  openDialog: PropTypes.bool.isRequired,
  exampleUploadData: PropTypes.object.isRequired
}

export default function ExampleUploadButton({...props}) {
  const {api} = useApi()
  const errors = useErrors()

  const [openDialog, setOpenDialog] = useState(false)

  const handleClickOpen = () => {
    setOpenDialog(true)
  }

  const handleClose = (value) => {
    (value && api.post(`/uploads?local_path=${value}`)
      .then(window.location.reload())
      .catch((error) => {
        errors.raiseError(error)
      })
    )
    setOpenDialog(false)
  }

  return (
    <div>
      <Button variant="contained" onClick={handleClickOpen} {...props}>
        Select from Example Uploads
      </Button>
      <ExampleUploadDialog
        exampleUploadData={exampleUploads}
        openDialog={openDialog}
        onClose={handleClose}
      />
    </div>
  )
}
ExampleUploadButton.propTypes = {
  isDisable: PropTypes.bool
}
