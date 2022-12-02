import React, { useCallback, useEffect, useMemo, useState } from 'react'
import {
  Box, Button, TextField, Dialog, DialogTitle, IconButton,
  DialogContent, DialogContentText, DialogActions, Tooltip
} from '@material-ui/core'
import { Autocomplete } from '@material-ui/lab'
import { useApi } from '../api'
import { useErrors } from '../errors'
import { getUrl } from '../nav/Routes'
import { useHistory, useLocation } from 'react-router-dom'
import { getUrlFromDefinition, SectionMDef, useGlobalMetainfo } from '../archive/metainfo'
import { useUploadPageContext } from './UploadPageContext'
import SearchIcon from '@material-ui/icons/Search'
import SectionSelectDialog from './SectionSelectDialog'

function addSearchIconToEndAdornment(endAdornment, onClick) {
  const children = React.Children.toArray(endAdornment.props.children)
  children.push(<IconButton size={'small'} onClick={onClick}>
    <Tooltip title="Search custom schema">
      <SearchIcon/>
    </Tooltip>
  </IconButton>)
  return React.cloneElement(endAdornment, {}, children)
}

const CreateEntry = React.memo(() => {
  const {deploymentUrl, uploadId, isProcessing} = useUploadPageContext()
  const {api} = useApi()
  const {raiseError} = useErrors()
  const globalMetainfo = useGlobalMetainfo()
  const [templates, setTemplates] = useState([])
  const [customTemplates, setCustomTemplates] = useState([])
  const [template, setTemplate] = useState()
  const [name, setName] = useState('')
  const [open, setOpen] = useState(false)
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
      const globalDefinitions = await globalMetainfo.getDefs()
      const globalTemplates = getTemplatesFromDefinitions(
        globalDefinitions, '__global__', null, section => section._qualifiedName)
      return globalTemplates.map(template => ({group: 'OASIS', ...template}))
    }

    getTemplates().then(setTemplates).catch(raiseError)
  }, [api, raiseError, setTemplates, globalMetainfo, isProcessing, deploymentUrl, uploadId])

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

  const handleChange = useCallback((event, value) => {
    setTemplate(value)
  }, [])

  const getTemplateFromDefinition = useCallback((dataSection, prefix, archive, getReference) => {
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
      entry_id: prefix,
      value: dataSection.name,
      archive: {
        data: {
          m_def: getReference(dataSection),
          ...template
        }
      }
    }
  }, [])

  const handleSelect = useCallback((value) => {
    const data = value.data
    let customTemplate = getTemplateFromDefinition(data.sectionDef, data.archive.metadata.entry_id, data.archive,
        section => {
          return getUrlFromDefinition(section, {deploymentUrl, uploadId}, true)
        })
    customTemplate = {group: 'Recent searches', ...customTemplate}
    setTemplate(customTemplate)
    if (!customTemplates.find(template => customTemplate?.id === template.id)) {
      setCustomTemplates([customTemplate, ...customTemplates])
    }
    setOpen(false)
  }, [customTemplates, getTemplateFromDefinition, deploymentUrl, uploadId])

  const allTemplates = useMemo(() => [...customTemplates, ...templates], [customTemplates, templates])
  const filtersLocked = useMemo(() => ({entry_type: ['Schema']}), [])

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
          onChange={event => setName(event.target.value.replace(/ /g, '_'))}/>
      </Box>
      <Autocomplete
        disabled={!allTemplates}
        value={allTemplates && allTemplates?.find(template => template.id === selectedTemplate?.id) ? selectedTemplate : null}
        onChange={handleChange}
        options={allTemplates || []}
        getOptionLabel={(option) => option.label}
        groupBy={(option) => option.group}
        style={{ width: 400 }}
        renderInput={(params) => {
          return (
            <TextField
              {...params}
              label="schema"
              variant="filled"
              InputProps={{
                ...params.InputProps,
                endAdornment: addSearchIconToEndAdornment(
                  params.InputProps.endAdornment,
                  () => setOpen(true)
                )
              }}
            />
          )
        }}
      />
      <SectionSelectDialog
          open={open}
          onCancel={() => setOpen(false)}
          onSelectedChanged={handleSelect}
          selected={customTemplates.find(customTemplate => customTemplate === template) && {entry_id: template?.entry_id, value: template?.value}}
          filtersLocked={filtersLocked}
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
