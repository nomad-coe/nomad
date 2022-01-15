import React, { useCallback, useContext, useState } from 'react'
import { Box, Button, TextField } from '@material-ui/core'
import { Autocomplete } from '@material-ui/lab'
import { uploadPageContext } from './UploadPage'
import { useApi } from '../api'
import { useErrors } from '../errors'

const templates = [
  {
    id: 'sample',
    label: 'Sample',
    description: 'A sample entry.',
    archive: {
      sample: [{
        name: 'unnamed'
      }]
    }
  },
  {
    id: 'experiment',
    label: 'Experiment',
    description: 'A sample experiment.',
    archive: {
      experiment: [{
        name: 'unnamed'
      }]
    }
  }
]

const CreateEntry = React.memo(function CreateEntry(props) {
  const {setUpload, data} = useContext(uploadPageContext)
  const uploadId = data.upload.upload_id
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [name, setName] = useState('')
  const [template, setTemplate] = useState(templates[0])

  const handleAdd = useCallback(() => {
    api.put(`uploads/${uploadId}/raw/?file_name=${name}.archive.json`, template.archive)
      .then(response => setUpload(response.data))
      .catch(raiseError)
  }, [api, raiseError, setUpload, template, name, uploadId])

  return <React.Fragment>
    <Box display="flex" flexDirection="row" alignItems="center" >
      <Box flexGrow={1} marginRight={2}>
        <TextField
          fullWidth variant="filled" label="name" value={name}
          onChange={event => setName(event.target.value.replaceAll(' ', '_'))}/>
      </Box>
      <Autocomplete
        value={template}
        onChange={(event, value) => setTemplate(value)}
        options={templates}
        getOptionLabel={(option) => option.label}
        style={{ width: 300 }}
        renderInput={(params) => <TextField {...params} label="template" variant="filled" />}
      />
    </Box>
    <Box display="flex" justifyContent="end" marginY={1}>
      <Button
        variant="contained" color="primary"
        disabled={!name || name === '' || !template}
        onClick={handleAdd}
      >
        Add
      </Button>
    </Box>
  </React.Fragment>
})

export default CreateEntry
