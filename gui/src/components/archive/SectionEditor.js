import React, { useCallback, useState } from 'react'
import PropTypes from 'prop-types'
import {Box, Button, makeStyles, Select, TextField, Typography} from '@material-ui/core'
import { useEntryContext } from '../entry/EntryContext'
import { useApi } from '../api'
import { useErrors } from '../errors'

const PropertyEditor = React.memo(function PropertyEditor({property, value, onChange}) {
  const handleChange = useCallback((value) => {
    if (onChange) {
      onChange(value)
    }
  }, [onChange])
  if (property.type?.type_kind === 'python' && property.type?.type_data === 'str' && property.shape?.length === 0) {
    return <TextField
      fullWidth variant="filled" value={value || ''} label={property.name}
      onChange={event => handleChange(event.target.value)}/>
  } else if (property.type?.type_kind === 'Enum' && property.shape?.length === 0) {
    return <Select variant="filled" native style={{width: '100%'}} inputProps={{name: property.name, id: property.name}}>
      {property.type?.type_data.map(item => <option key={item}>{item}</option>)}
    </Select>
  }
  return ''
})
PropertyEditor.propTypes = {
  property: PropTypes.object.isRequired,
  value: PropTypes.any,
  onChange: PropTypes.func
}

const useSectionEditorStyles = makeStyles(theme => ({
  root: {
    minWidth: 600
  }
}))
const SectionEditor = React.memo(function SectionEditor({sectionDef, section, onChange}) {
  const classes = useSectionEditorStyles()
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [saving, setSaving] = useState(false)
  const [version, setVersion] = useState(0)
  const [savedVersion, setSavedVersion] = useState(0)
  const {metadata, archive, reload} = useEntryContext()
  const hasChanges = version > savedVersion
  const handleChange = useCallback((property, value) => {
    section[property.name] = value
    if (onChange) {
      onChange(section)
    }
    setVersion(value => value + 1)
  }, [section, onChange, setVersion])
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
          value={section && section[property.name]} onChange={value => handleChange(property, value)}
        />
      </Box>
    ))}
    <Box display="flex" flexDirection="row" alignItems="center" marginY={1}>
      <Box flexGrow={1}>
        <Typography>{hasChanges ? 'Not yet saved' : 'All changes saved'}</Typography>
      </Box>
      <Button
        color="primary" variant="contained"
        disabled={!hasChanges || saving} onClick={handleSave}
      >
        {saving ? 'Saving...' : 'Save'}
      </Button>
    </Box>
  </div>
})
SectionEditor.propTypes = {
  sectionDef: PropTypes.object.isRequired,
  section: PropTypes.object,
  onChange: PropTypes.func,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
}

export default SectionEditor
