import React, { useCallback, useEffect, useState } from 'react'
import { Box, Button, TextField } from '@material-ui/core'
import { Autocomplete } from '@material-ui/lab'
import { useApi } from '../api'
import { useErrors } from '../errors'
import { getUrl } from '../nav/Routes'
import { useHistory, useLocation } from 'react-router-dom'
import { createMetainfo, getSectionReference, SectionMDef, useGlobalMetainfo } from '../archive/metainfo'
import { useUploadContext } from './UploadContext'

const CreateEntry = React.memo(function CreateEntry(props) {
  const {data} = useUploadContext()
  const {api} = useApi()
  const {raiseError} = useErrors()
  const globalMetainfo = useGlobalMetainfo()
  const [templates, setTemplates] = useState()
  const [template, setTemplate] = useState()
  const uploadId = data.upload.upload_id
  const [name, setName] = useState('')
  const history = useHistory()
  const location = useLocation()
  const selectedTemplate = template || templates?.[0] || null

  useEffect(() => {
    // TODO this whole thing is quite expensive to repeat on each upload page?
    // There needs to be a register/cache based on hashes or something
    if (!globalMetainfo) {
      return
    }

    const getTemplatesFromDefinitions = (definitions, prefix, getReference) => {
      return definitions.filter(definition => {
        if (definition.m_def !== SectionMDef) {
          return false
        }
        return definition.base_sections?.find(baseSection => baseSection.name === 'EntryData')
      }).map(dataSection => ({
        id: `${prefix}:${dataSection._qualifiedName}`,
        label: dataSection.name,
        archive: {
          data: {
            m_def: getReference(dataSection)
          }
        }
      }))
    }

    const getTemplates = async () => {
      // TODO paginate?
      const response = await api.post(`entries/archive/query`, {
        owner: 'visible',
        query: {
          quantities: 'definitions.section_definitions'
        },
        required: {
          definitions: '*',
          metadata: {
            mainfile: '*'
          }
        }
      })
      for (const data of response.data) {
        const archive = data.archive
        // TODO we should not createMetainfo all the time? There needs to be a
        // register/cache based on hashes or something
        archive._metainfo = await createMetainfo(archive, globalMetainfo, {api: api, archive: archive})
      }
      const customTemplates = response.data.reduce((templates, data) => {
        const archive = data.archive
        const newTemplates = getTemplatesFromDefinitions(
          archive.definitions.section_definitions, data.entry_id,
          section => {
            const fragment = getSectionReference(section)
            return `../uploads/${data.upload_id}/raw/${archive.metadata.mainfile}#/definitions${fragment}`
          })
        newTemplates.forEach(template => {
          templates.push(template)
        })
        return templates
      }, [])
      const globalDefinitions = await globalMetainfo.getDefs()
      const globalTemplates = getTemplatesFromDefinitions(
        globalDefinitions, '__global__', section => section._qualifiedName)
      return customTemplates.concat(globalTemplates)
    }

    getTemplates().then(setTemplates).catch(raiseError)
  }, [api, raiseError, setTemplates, globalMetainfo])

  const handleAdd = useCallback(() => {
    api.put(`uploads/${uploadId}/raw/?file_name=${name}.archive.json&wait_for_processing=true`, selectedTemplate.archive)
      .then(response => {
        // TODO handle processing errors
        const entryId = response.processing.entry_id
        history.push(getUrl(`entry/id/${entryId}`, location))
      })
      .catch(raiseError)
  }, [api, raiseError, selectedTemplate, name, uploadId, history, location])

  return <React.Fragment>
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
        style={{ width: 300 }}
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
