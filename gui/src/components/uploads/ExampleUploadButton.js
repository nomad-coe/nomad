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
import { Button, Dialog, DialogTitle } from '@material-ui/core'
import React, { useState } from 'react'
import PropTypes from 'prop-types'
import { makeStyles } from '@material-ui/core/styles'
import List from '@material-ui/core/List'
import ListItem from '@material-ui/core/ListItem'
import ListItemText from '@material-ui/core/ListItemText'
import { blue } from '@material-ui/core/colors'
import CloudUploadIcon from '@material-ui/icons/CloudUpload'
import { useHistory, useLocation } from 'react-router-dom'
import { useApi } from '../api'
import { useErrors } from '../errors'

const exampleUploads = {
  Theory: [
    {
      title: 'VASP Sample Uploads',
      description: 'Here goes the VASP description Here goes the VASP description Here goes the VASP description Here goes the VASP description Here goes the VASP description ',
      filePath: 'tests/data/proc/examples_vasp.zip'
    },
    {
      title: 'Smaller VASP Sample Uploads',
      description: 'Here goes the Smaller VASP description',
      filePath: 'tests/data/proc/examples_vasp.zip'
    }
  ],
  Tabular: [
    {
      title: 'ELN Sample Uploads',
      description: 'Here goes the ELN description',
      filePath: 'tests/data/proc/examples_vasp.zip'
    }
  ]
}

const useStyles = makeStyles({
  buttons: {
    backgroundColor: blue[100],
    color: blue[600]
  }
})

const ExampleUploadDialog = React.memo(function ExampleUploadDialog(props) {
  const { exampleUploadData, onClose, selectedUpload, open } = props
  const classes = useStyles()

  const handleClose = () => {
    onClose(selectedUpload)
  }

  return (
    <Dialog
      onClose={handleClose}
      aria-labelledby="Uploads-dialogBox"
      open={open}
      fullWidth
    >
      <DialogTitle id="Uploads-dialogBox">Select a sample Upload</DialogTitle>
      <List>
        {Object.keys(exampleUploadData).map((exampleType) => exampleUploadData[exampleType].map((uploadEntry, i) =>
          <ListItem key={i}>
            <ListItemText
              primary={`${uploadEntry.title}  (${exampleType})`}
              secondary={uploadEntry.description}
              style={{maxWidth: '450px'}}
            />
            <Button
              variant="contained"
              color="primary"
              className={classes.button}
              startIcon={<CloudUploadIcon />}
              onClick={() => onClose(uploadEntry.filePath)}
              style={{marginLeft: '20px'}}
            >
              Select
            </Button>
          </ListItem>)
        )}
      </List>
    </Dialog>
  )
})

ExampleUploadDialog.propTypes = {
  onClose: PropTypes.func.isRequired,
  open: PropTypes.bool.isRequired,
  exampleUploadData: PropTypes.object.isRequired,
  selectedUpload: PropTypes.string
}

export default function ExampleUploadButton({...props}) {
  const {api} = useApi()
  const errors = useErrors()
  const history = useHistory()
  const location = useLocation()

  const [selectedUpload, setSelectedUpload] = useState(null)
  const [open, setOpen] = useState(false)

  const handleClickOpen = () => {
    setOpen(true)
  }

  const handleClose = (value) => {
    api.post(`/uploads?local_path=${value}`)
      .then(setSelectedUpload(value)
        // (upload) => {
        //   history.push(getUrl(`upload/id/${upload.upload_id}`, location))
        // }
      )
      .catch((error) => {
        errors.raiseError(error)
      })
    // setSelectedUpload(value)
    setOpen(false)
    console.log(selectedUpload)
  }

  return (
    <div>
      <Button variant="contained" onClick={handleClickOpen} {...props}>
        Select from Example Uploads
      </Button>
      <ExampleUploadDialog
        exampleUploadData={exampleUploads}
        selectedUpload={selectedUpload}
        open={open}
        onClose={handleClose}
        api={api}
        location={location}
        history={history}
        errors={errors}
      />
    </div>
  )
}
ExampleUploadButton.propTypes = {
  isDisable: PropTypes.bool
}
