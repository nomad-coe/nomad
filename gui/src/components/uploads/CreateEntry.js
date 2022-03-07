import React, { useCallback, useMemo, useState } from 'react'
import { Box, Button, TextField } from '@material-ui/core'
import { Autocomplete } from '@material-ui/lab'
import { useApi } from '../api'
import { useErrors } from '../errors'
import { getUrl } from '../nav/Routes'
import { useHistory, useLocation } from 'react-router-dom'
import { defs as metainfoDefs, SectionMDef, resolveRef } from '../archive/metainfo'
import { useUploadContext } from './UploadContext'

const CreateEntry = React.memo(function CreateEntry(props) {
  const {data} = useUploadContext()
  const templates = useMemo(() => {
    return metainfoDefs
      .filter(definition => {
        if (definition.m_def !== SectionMDef) {
          return false
        }
        return definition.base_sections?.find(baseSection => resolveRef(baseSection).name === 'EntryData')
      })
      .map(dataSection => ({
        id: dataSection._qualifiedName,
        label: dataSection.name,
        archive: {
          data: {
            m_def: dataSection._qualifiedName
          }
        }
      }))
  }, [])
  const uploadId = data.upload.upload_id
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [name, setName] = useState('')
  const [template, setTemplate] = useState(templates[0])
  const history = useHistory()
  const location = useLocation()

  const handleAdd = useCallback(() => {
    api.put(`uploads/${uploadId}/raw/?file_name=${name}.archive.json&wait_for_processing=true`, template.archive)
      .then(response => {
        // TODO handle processing errors
        const entryId = response.processing.entry_id
        history.push(getUrl(`entry/id/${entryId}`, location))
      })
      .catch(raiseError)
  }, [api, raiseError, template, name, uploadId, history, location])

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
        renderInput={(params) => <TextField {...params} label="type" variant="filled" />}
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
