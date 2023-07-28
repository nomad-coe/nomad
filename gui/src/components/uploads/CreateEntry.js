import React, { useCallback, useEffect, useMemo, useState } from 'react'
import {
  Button, TextField, Dialog, DialogTitle, IconButton, DialogContent,
  DialogContentText, DialogActions, Tooltip, makeStyles, Box, RadioGroup, FormControlLabel, Radio,
  FormControl
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
import SectionSelectAutocomplete from './SectionSelectAutocomplete'

const useStyles = makeStyles(theme => ({
  dialog: {
    width: '100%',
    minWidth: 600
  },
  nameField: {
    marginBottom: theme.spacing(3),
    width: '100%'
  },
  schemaField: {
    marginBottom: theme.spacing(1),
    width: '100%'
  },
  radioGroup: {
    marginTop: 0,
    display: 'flex'
  }
}))

function addSearchIconToEndAdornment(endAdornment, onClick) {
  const children = React.Children.toArray(endAdornment.props.children)
  children.push(<IconButton key={'searchButton'} size={'small'} onClick={onClick}>
    <Tooltip title="Search custom schema">
      <SearchIcon/>
    </Tooltip>
  </IconButton>)
  return React.cloneElement(endAdornment, {}, children)
}

const CreateEntry = React.memo((props) => {
  const classes = useStyles()
  const {deploymentUrl, uploadId, isProcessing} = useUploadPageContext()
  const {api} = useApi()
  const {raiseError} = useErrors()
  const globalMetainfo = useGlobalMetainfo()
  const [builtInTemplates, setBuiltInTemplates] = useState([])
  const [builtInTemplate, setBuiltInTemplate] = useState(null)
  const [customTemplate, setCustomTemplate] = useState(null)
  const [name, setName] = useState('')
  const [open, setOpen] = useState(false)
  const history = useHistory()
  const location = useLocation()
  const [openEntryAlreadyExistsDialog, setOpenEntryAlreadyExistsDialog] = useState(false)
  const [openCreateEntryDialog, setOpenCreateEntryDialog] = useState(false)
  const [schemaType, setSchemaType] = useState('built-in')
  const [error, setError] = useState()

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

    const getTemplatesFromGlobalDefinitions = (definitions, prefix, archive, getReference) => {
      const templates = definitions.filter(definition => {
        if (definition.m_def !== SectionMDef) {
          return false
        }
        return definition._allBaseSections?.find(baseSection => baseSection.name === 'EntryData')
      }).map(definition => {
        const label = definition.label || definition.name
        const template = definition.m_annotations?.template?.[0] || {}
        const findCategory = definition => (
          definition.categories.find(category => (
            category.categories.find(category => category._qualifiedName === 'nomad.datamodel.data.EntryDataCategory')
          ))
        )
        let category
        for (const baseDefinition of [definition, ...definition._allBaseSections]) {
          category = findCategory(baseDefinition)
          if (category) {
            break
          }
        }
        return {
          id: `${prefix}:${definition._qualifiedName}`,
          categoryLabel: category?.label || category?.name || 'Other',
          categoryName: category?.name || 'Other',
          definition,
          label,
          archive: {
            data: {
              m_def: getReference(definition),
              ...template
            }
          }
        }
      })
      const categoryNames = templates.map(template => template.categoryName)
      const categoryOrder = ['BasicElnCategory', 'ElnExampleCategory', 'ElnIntegrationCategory', 'UseCaseElnCategory', 'WorkflowsElnCategory', 'Other']
      categoryNames.sort((a, b) => categoryOrder.indexOf(a) - categoryOrder.indexOf(b))
      templates.sort((a, b) => a.label.localeCompare(b.label))
      templates.sort((a, b) => categoryNames.indexOf(a.categoryName) - categoryNames.indexOf(b.categoryName))
      return templates
    }

    const getTemplates = async () => {
      const globalDefinitions = await globalMetainfo.getDefs()
      const globalTemplates = getTemplatesFromGlobalDefinitions(
        globalDefinitions, '__global__', null, section => section._qualifiedName)
      return globalTemplates.map(template => ({group: 'OASIS', ...template}))
    }

    getTemplates().then(templates => {
      setBuiltInTemplates(templates)
      setBuiltInTemplate(templates[0])
    }).catch(raiseError)
  }, [api, raiseError, setBuiltInTemplates, globalMetainfo, isProcessing, deploymentUrl, uploadId])

  const handleAdd = useCallback(() => {
    const selectedTemplate = schemaType === 'built-in' ? builtInTemplate : customTemplate
    api.put(`uploads/${uploadId}/raw/?file_name=${name}.archive.json&overwrite_if_exists=false&wait_for_processing=true`, selectedTemplate.archive)
      .then(response => {
        // TODO handle processing errors
        const entryId = response.processing.entry_id
        history.push(getUrl(`entry/id/${entryId}/data/data`, location))
        setOpenCreateEntryDialog(false)
      })
      .catch(error => {
        if (error.apiMessage?.startsWith('The provided path already exists')) {
          setOpenEntryAlreadyExistsDialog(true)
        } else {
          raiseError(error)
        }
      })
  }, [setOpenEntryAlreadyExistsDialog, api, raiseError, schemaType, builtInTemplate, customTemplate, name, uploadId, history, location])

  const handleBuiltInChange = useCallback((event, value) => {
    setBuiltInTemplate(value)
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

  const handleCustomChange = useCallback((value) => {
    const data = value?.data
    if (!data) {
      setCustomTemplate(null)
      return
    }
    const customTemplate = getTemplateFromDefinition(data.sectionDef, data.archive.metadata.entry_id, data.archive,
      section => {
        return getUrlFromDefinition(section, {deploymentUrl, uploadId}, true)
      })
    setCustomTemplate(customTemplate)
  }, [deploymentUrl, getTemplateFromDefinition, uploadId])

  const handleSelect = useCallback((value) => {
    const data = value?.data
    if (!data) {
      setCustomTemplate(null)
      return
    }
    const customTemplate = getTemplateFromDefinition(data.sectionDef, data.archive.metadata.entry_id, data.archive,
        section => {
          return getUrlFromDefinition(section, {deploymentUrl, uploadId}, true)
        })
    setCustomTemplate(customTemplate)
    setOpen(false)
  }, [getTemplateFromDefinition, deploymentUrl, uploadId])

  const filtersLocked = useMemo(() => ({'section_defs.definition_qualified_name': ['nomad.metainfo.metainfo.Definition']}), [])

  const handleChangeTab = (event) => {
    setSchemaType(event.target.value)
  }

  const handleError = useCallback((error) => {
    setError(error)
  }, [])

  return <React.Fragment>
    <Button
      onClick={() => setOpenCreateEntryDialog(true)}
      disabled={isProcessing}
      {...props}
    />
    <Dialog
      classes={{paper: classes.dialog}}
      open={openCreateEntryDialog}
      onClose={() => setOpenCreateEntryDialog(false)}
      data-testid='create-entry-dialog'
    >
      <DialogTitle>Create new entry from schema</DialogTitle>
      <DialogContent>
        <DialogContentText>
          Name of the entry
        </DialogContentText>
        <TextField
          className={classes.nameField}
          variant="filled"
          label="name"
          value={name}
          onChange={event => setName(event.target.value.replace(/ /g, '_'))}/>
        <DialogContentText>
          Select a schema
        </DialogContentText>
        <FormControl component="fieldset" className={classes.radioGroup}>
          <RadioGroup
            row
            aria-label="position"
            name="position"
            value={schemaType}
            onChange={handleChangeTab}
            className={classes.radioGroup}
          >
            <FormControlLabel
              value="built-in"
              control={<Radio color="primary" />}
              label="Built-in schema"
              labelPlacement="end"
            />
            <FormControlLabel
              value="custom"
              control={<Radio color="primary" data-testid='custom-schema-radio'/>}
              label="Custom schema"
              labelPlacement="end"
            />
          </RadioGroup>
        </FormControl>
        {schemaType === 'built-in' && <Autocomplete
          className={classes.schemaField}
          disabled={!builtInTemplates}
          value={builtInTemplates && builtInTemplates?.find(template => template.id === builtInTemplate?.id) ? builtInTemplate : null}
          onChange={handleBuiltInChange}
          options={builtInTemplates || []}
          groupBy={option => option.categoryLabel}
          getOptionLabel={option => option.label}
          renderInput={(params) => {
            return (
                <TextField
                  {...params}
                  label='built-in schema'
                  variant='filled'
                  data-testid='builtin-select-schema'
                />
              )
          }}
        />}
        {schemaType === 'custom' && <Box>
          <SectionSelectAutocomplete
              onValueChanged={handleCustomChange}
              value={customTemplate}
              filtersLocked={filtersLocked}
              onError={handleError}
              renderInput={(params) => {
                return (
                    <TextField
                        {...params}
                        variant='filled'
                        error={!!error}
                        helperText={error}
                        label='custom schema'
                        placeholder={'search by entry name or file name'}
                        InputProps={{
                          ...params.InputProps,
                          endAdornment: addSearchIconToEndAdornment(
                              params.InputProps.endAdornment,
                              () => setOpen(true)
                          )
                        }}
                        data-testid='custom-select-schema'
                    />
                )
              }}
          />
        </Box>}
        <DialogContentText>
          {name ? `File name: ${name}.archive.json` : ''}
        </DialogContentText>
        <SectionSelectDialog
          open={open}
          onCancel={() => setOpen(false)}
          onSelectedChanged={handleSelect}
          selected={customTemplate && {entry_id: customTemplate?.entry_id, value: customTemplate?.value}}
          filtersLocked={filtersLocked}
        />
      </DialogContent>
      <DialogActions>
        <span style={{flexGrow: 1}} />
        <Button onClick={() => setOpenCreateEntryDialog(false)} color="secondary">
          Cancel
        </Button>
        <Button disabled={!name || name === '' || !(schemaType === 'built-in' ? builtInTemplate : customTemplate)} onClick={handleAdd} color="secondary">
          Create
        </Button>
      </DialogActions>
    </Dialog>
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
  </React.Fragment>
})

export default CreateEntry
