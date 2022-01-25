import React, { useCallback, useState } from 'react'
import PropTypes from 'prop-types'
import {Box, Button, FormControl, InputLabel, makeStyles, Select, TextField, Typography, InputAdornment} from '@material-ui/core'
import { useEntryContext } from '../entry/EntryContext'
import { useApi } from '../api'
import { useErrors } from '../errors'
import {Compartment} from './Browser'
import _ from 'lodash'

const useStyles = makeStyles(theme => ({
  adornment: {
    marginRight: theme.spacing(3)
  }
}))

const PropertyEditor = React.memo(function PropertyEditor({property, section, value, nestedPath, parentVersion, setParentVersion, onChange}) {
  const classes = useStyles()
  const [validationError, setValidationError] = useState('')
  const handleChange = useCallback((value) => {
    if (onChange) {
      onChange(value)
    }
  }, [onChange])

  function isFloat(str) {
    const num = Number(str)
    return !!num
  }

  function isInteger(str) {
    const num = Number(str)
    return Number.isInteger(num)
  }

  function isUInteger(str) {
    const num = Number(str)
    return Number.isInteger(num) && num > 0
  }

  const handleFloatValidator = (event) => {
    (isFloat(event.target.value) || event.target.value === '' ? setValidationError('') : setValidationError('Please enter a valid number!'))
  }

  const handleIntegerValidator = (event) => {
    (isInteger(event.target.value) || event.target.value === '' ? setValidationError('') : setValidationError('Please enter an integer number!'))
  }

  const handleUIntegerValidator = (event) => {
    (isUInteger(event.target.value) || event.target.value === '' ? setValidationError('') : setValidationError('Please enter an unsigned integer number!'))
  }

  if (property.m_def === 'SubSection') {
    let currentPath = (nestedPath ? `${nestedPath}.${property.name}` : property.name)
    return <Compartment title={currentPath}>
      <SectionEditor sectionDef={property._subSection} section={section} nestedPath={currentPath} parentVersion={parentVersion} setParentVersion={setParentVersion}/>
    </Compartment>
  } else if (property.type?.type_kind === 'python') {
    if (property.type?.type_data === 'str' && property.shape?.length === 0) {
      return <TextField
        fullWidth variant="filled" size='small' value={value || ''} label={property.name} multiline={property.name === 'description'} minRows={4}
        InputProps={{endAdornment: <InputAdornment className={classes.adornment} position='end'>{property.unit}</InputAdornment>}}
        onChange={event => handleChange(event.target.value)}/>
    }
  } else if (property.type?.type_kind === 'numpy') {
    if (property.type?.type_data === 'int64' && property.shape?.length === 0) {
      return <TextField onBlur={handleIntegerValidator} error={!!validationError} helperText={validationError}
        fullWidth variant="filled" size='small' value={value || ''} label={property.name}
        InputProps={{endAdornment: <InputAdornment className={classes.adornment} position='end'>{property.unit}</InputAdornment>}}
        onChange={event => handleChange(event.target.value)}/>
    } else if (property.type?.type_data === 'uint64' && property.shape?.length === 0) {
      return <TextField onBlur={handleUIntegerValidator} error={!!validationError} helperText={validationError}
        fullWidth variant="filled" size='small' value={value || ''} label={property.name}
        InputProps={{endAdornment: <InputAdornment className={classes.adornment} position='end'>{property.unit}</InputAdornment>}}
        onChange={event => handleChange(event.target.value)}/>
    } else if (property.type?.type_data === 'float64' && property.shape?.length === 0) {
      return <TextField onBlur={handleFloatValidator} error={!!validationError} helperText={validationError}
        fullWidth variant="filled" size='small' value={value || ''} label={property.name}
        InputProps={{endAdornment: <InputAdornment className={classes.adornment} position='end'>{property.unit}</InputAdornment>}}
        onChange={event => handleChange(event.target.value)}/>
    }
  } else if (property.type?.type_kind === 'Enum' && property.shape?.length === 0) {
    return <FormControl variant='filled' size='small' fullWidth>
      <InputLabel htmlFor={property.name} shrink={value !== undefined && value !== ''}>{property.name}</InputLabel>
      <Select native value={value}
        endAdornment={<InputAdornment className={classes.adornment} position='end'>{property.unit}</InputAdornment>}
        onChange={(event) => handleChange(event.target.value)}>
        <option key={''}>{''}</option>
        {property.type?.type_data.map(item => <option key={item}>{item}</option>)}
      </Select>
    </FormControl>
  }
  return ''
})
PropertyEditor.propTypes = {
  property: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  value: PropTypes.any,
  parentVersion: PropTypes.bool,
  setParentVersion: PropTypes.func,
  nestedPath: PropTypes.string,
  onChange: PropTypes.func
}

const useSectionEditorStyles = makeStyles(theme => ({
  root: {
    minWidth: 600
  }
}))
const SectionEditor = React.memo(function SectionEditor({sectionDef, section, nestedPath, parentVersion, setParentVersion, onChange}) {
  const classes = useSectionEditorStyles()
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [saving, setSaving] = useState(false)
  const [version, setVersion] = useState(0)
  const [savedVersion, setSavedVersion] = useState(0)
  const {metadata, archive, reload} = useEntryContext()
  const hasChanges = version > savedVersion

  const handleChange = useCallback((property, value) => {
    let currentPath = (nestedPath ? `${nestedPath}.${property.name}` : property.name)
    _.set(section, currentPath, value)
    if (onChange) {
      onChange(section)
    }
    if (parentVersion === undefined) {
      setVersion(value => value + 1)
    } else {
      setParentVersion(value => value + 1)
    }
  }, [section, onChange, setVersion, nestedPath, parentVersion, setParentVersion])

  const handleSave = useCallback(() => {
    const uploadId = metadata.upload_id
    const {mainfile} = metadata
    if (uploadId) {
      const separatorIndex = mainfile.lastIndexOf('/')
      const path = mainfile.substring(0, separatorIndex + 1)
      const fileName = mainfile.substring(separatorIndex + 1)
      const newArchive = {...archive}
      delete newArchive.metadata
      delete newArchive.results
      delete newArchive.processing_logs
      api.put(`/uploads/${uploadId}/raw/${path}?file_name=${fileName}`, newArchive)
        .then(() => {
          // TODO this is just a hack, wait a bit for reprocessing too complete
          window.setTimeout(reload, 500)
        })
        .catch(raiseError)
        .finally(() => setSaving(false))
      setSaving(true)
    }
    setSavedVersion(version)
  }, [setSavedVersion, version, api, raiseError, setSaving, archive, metadata, reload])

  return <div className={classes.root}>
    {sectionDef._allProperties.map(property => (
      <Box marginBottom={1} key={property.name}>
        <PropertyEditor
          property={property}
          section={section}
          value={section && _.get(section, (nestedPath ? `${nestedPath}.${property.name}` : property.name))} onChange={value => handleChange(property, value)}
          nestedPath={nestedPath}
          parentVersion={(parentVersion || version)}
          setParentVersion={(setParentVersion || setVersion)}
        />
      </Box>
    ))}
    {parentVersion === undefined && <Box display="flex" flexDirection="row" alignItems="center" marginY={1}>
      <Box flexGrow={1}>
        <Typography>{hasChanges ? 'Not yet saved' : 'All changes saved'}</Typography>
      </Box>
      <Button
        color="primary" variant="contained"
        disabled={!hasChanges || saving} onClick={handleSave}
      >
        {saving ? 'Saving...' : 'Save'}
      </Button>
    </Box>}
  </div>
})
SectionEditor.propTypes = {
  sectionDef: PropTypes.object.isRequired,
  section: PropTypes.object,
  nestedPath: PropTypes.string,
  onChange: PropTypes.func,
  parentVersion: PropTypes.bool,
  setParentVersion: PropTypes.func,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
}

export default SectionEditor
