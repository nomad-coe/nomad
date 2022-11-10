import React, { useCallback, useEffect, useState } from 'react'
import { Box, Button, TextField, Dialog, DialogTitle, DialogContent, DialogContentText, DialogActions } from '@material-ui/core'
import { Autocomplete } from '@material-ui/lab'
import { useApi } from '../api'
import { useErrors } from '../errors'
import { getUrl } from '../nav/Routes'
import { useHistory, useLocation } from 'react-router-dom'
import { getUrlFromDefinition, SectionMDef, useGlobalMetainfo } from '../archive/metainfo'
import { useUploadPageContext } from './UploadPageContext'

const CreateEntry = React.memo(function CreateEntry(props) {
  const {installationUrl, uploadId, isProcessing} = useUploadPageContext()
  const {api} = useApi()
  const {raiseError} = useErrors()
  const globalMetainfo = useGlobalMetainfo()
  const [templates, setTemplates] = useState([])
  const [template, setTemplate] = useState()
  const [name, setName] = useState('')
  const history = useHistory()
  const location = useLocation()
  const selectedTemplate = template || templates?.[0] || null
  const [openEntryAlreadyExistsDialog, setOpenEntryAlreadyExistsDialog] = useState(false)

  useEffect(() => {
    // TODO this whole thing is quite expensive to repeat on each upload page?
    // There needs to be a register/cache based on hashes or something
    if (!globalMetainfo) {
      return
    }

    // Do not reload all possible entries while still processing.
    if (isProcessing) {
      return
    }

    const getTemplatesFromDefinitions = (definitions, prefix, archive, getReference) => {
      return definitions.filter(definition => {
        if (definition.m_def !== SectionMDef) {
          return false
        }
        return definition._allBaseSections?.find(baseSection => baseSection.name === 'EntryData')
      }).map(dataSection => {
        let label = dataSection.name
        const entry_name = archive?.metadata?.entry_name
        const mainfile = archive?.metadata?.mainfile

        if (entry_name) {
          label = `${label} (${entry_name})`
        } else if (mainfile) {
          const file_name = mainfile.substring(mainfile.lastIndexOf('/') + 1)
          label = `${label} (${file_name})`
        } else if (prefix === '__global__') {
          label = `${label} (OASIS)`
        }
        const template = dataSection.m_annotations?.template?.[0] || {}
        return {
          id: `${prefix}:${dataSection._qualifiedName}`,
          label: label,
          archive: {
            data: {
              m_def: getReference(dataSection),
              ...template
            }
          }
        }
      })
    }

    const getTemplates = async () => {
      const customMetainfos = await globalMetainfo.fetchAllCustomMetainfos(true)
      const customTemplates = customMetainfos.reduce((templates, metainfo) => {
        const archive = metainfo._data
        const newTemplates = getTemplatesFromDefinitions(
          archive.definitions.section_definitions, archive.metadata.entry_id, archive,
          section => {
            return getUrlFromDefinition(section, {installationUrl, uploadId}, true)
          })
        newTemplates.forEach(template => {
          templates.push(template)
        })
        return templates
      }, [])
      const globalDefinitions = await globalMetainfo.getDefs()
      const globalTemplates = getTemplatesFromDefinitions(
        globalDefinitions, '__global__', null, section => section._qualifiedName)
      return customTemplates.concat(globalTemplates)
    }

    getTemplates().then(setTemplates).catch(raiseError)
  }, [api, raiseError, setTemplates, globalMetainfo, isProcessing, installationUrl, uploadId])

  const handleAdd = useCallback(() => {
    api.put(`uploads/${uploadId}/raw/?file_name=${name}.archive.json&overwrite_if_exists=false&wait_for_processing=true`, selectedTemplate.archive)
      .then(response => {
        // TODO handle processing errors
        const entryId = response.processing.entry_id
        history.push(getUrl(`entry/id/${entryId}/data/data`, location))
      })
      .catch(error => {
        if (error.apiMessage?.startsWith('The provided path already exists')) {
          setOpenEntryAlreadyExistsDialog(true)
        } else {
          raiseError(error)
        }
      })
  }, [setOpenEntryAlreadyExistsDialog, api, raiseError, selectedTemplate, name, uploadId, history, location])

  return <React.Fragment>
    <Dialog
      open={openEntryAlreadyExistsDialog}
      onClose={() => setOpenEntryAlreadyExistsDialog(false)}
    >
      <DialogTitle>Entry already exists</DialogTitle>
      <DialogContent>
        <DialogContentText>
          A mainfile with the specified name already exists. Please choose a different name.
        </DialogContentText>
      </DialogContent>
      <DialogActions>
        <Button onClick={() => setOpenEntryAlreadyExistsDialog(false)} autoFocus>OK</Button>
      </DialogActions>
    </Dialog>
    <Box display="flex" flexDirection="row" alignItems="center" >
      <Box flexGrow={1} marginRight={2}>
        <TextField
          fullWidth variant="filled" label="name" value={name}
          onChange={event => setName(event.target.value.replaceAll(' ', '_'))}/>
      </Box>
      <Autocomplete
        disabled={!templates}
        value={selectedTemplate}
        onChange={(event, value) => setTemplate(value)}
        options={templates || []}
        getOptionLabel={(option) => option.label}
        style={{ width: 400 }}
        renderInput={(params) => <TextField {...params} label="type" variant="filled" />}
      />
    </Box>
    <Box display="flex" justifyContent="end" marginY={1}>
      <Button
        variant="contained" color="primary"
        disabled={!name || name === '' || !selectedTemplate}
        onClick={handleAdd}
      >
        Add
      </Button>
    </Box>
  </React.Fragment>
})

export default CreateEntry
