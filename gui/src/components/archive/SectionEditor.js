import React, {useCallback, useState, useContext} from 'react'
import PropTypes from 'prop-types'
import {
  Box,
  Button,
  FormControl,
  InputLabel,
  makeStyles,
  Select,
  TextField,
  Typography,
  InputAdornment,
  DialogContent, Dialog
} from '@material-ui/core'
import { useEntryContext } from '../entry/EntryContext'
import { useApi } from '../api'
import { useErrors } from '../errors'
import _ from 'lodash'
import {formatSubSectionName, Item, laneContext} from './Browser'
import { SubSectionList } from './ArchiveBrowser'
import DialogContentText from '@material-ui/core/DialogContentText'
import DialogActions from '@material-ui/core/DialogActions'

const useStyles = makeStyles(theme => ({
  adornment: {
    marginRight: theme.spacing(3)
  }
}))

const PropertyEditor = React.memo(function PropertyEditor({property, section, value, nestedPath, onChange}) {
  const lane = useContext(laneContext)
  const classes = useStyles()
  const [validationError, setValidationError] = useState('')
  const handleChange = useCallback((value) => {
    if (onChange) {
      onChange(value)
    }
  }, [onChange])

  const handleCreate = useCallback(() => {
    if (onChange) {
      onChange((property.repeats ? [{}] : {}))
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

  if (!lane) return ''

  if (property.m_def === 'SubSection') {
    // let currentPath = (nestedPath ? `${nestedPath}.${property.name}` : property.name)
    const key = property.name
    const isEmpty = section && section[key] === undefined
    if (isEmpty) {
      return <Item key={key} itemKey={key} disabled={true}>
        <Typography component="span">
          <Box fontWeight="bold" component="span">
            {formatSubSectionName(property.name)}
            <Button color='primary' size='small' onClick={() => handleCreate()}>
              Create
            </Button>
          </Box>
        </Typography>
      </Item>
    }
    if (property.repeats) {
      return <SubSectionList
        key={property.name}
        subSectionDef={property}/>
      // return <Item key={key} itemKey={key}>
      //   <Typography component="span">
      //     <Box fontWeight="bold" component="span">
      //       {formatSubSectionName(property.name)}
      //     </Box>
      //   </Typography>
      // </Item>
    } else {
      return <Item key={key} itemKey={key}>
        <Typography component="span">
          <Box fontWeight="bold" component="span">
            {formatSubSectionName(property.name)}
          </Box>
        </Typography>
      </Item>
    }
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
  nestedPath: PropTypes.string,
  onChange: PropTypes.func
}

const useSectionEditorStyles = makeStyles(theme => ({
  root: {
    minWidth: '400px',
    width: '100%'
  }
}))
const SectionEditor = React.memo(function SectionEditor({sectionDef, section, parent, nestedPath, onChange}) {
  const classes = useSectionEditorStyles()
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [saving, setSaving] = useState(false)
  const [version, setVersion] = useState(0)
  const [savedVersion, setSavedVersion] = useState(0)
  const {metadata, archive, reload} = useEntryContext()
  const [openConfirmDialog, setOpenConfirmDialog] = useState(false)
  const hasChanges = version > savedVersion

  const handleChange = useCallback((property, value) => {
    let currentPath = (nestedPath ? `${nestedPath}.${property.name}` : property.name)
    _.set(section, currentPath, value)
    if (onChange) {
      onChange(section)
    }
    setVersion(value => value + 1)
  }, [section, onChange, setVersion, nestedPath])

  const handleConfirm = () => {
    setOpenConfirmDialog(true)
  }

  const handleDelete = () => {
    section = parent
    _.set(section, sectionDef.name.toLowerCase(), undefined)
    if (onChange) {
      onChange(section)
    }
    setOpenConfirmDialog(false)
    handleSave()
  }

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
        />
      </Box>
    ))}
    <Box display="flex" flexDirection="row" alignItems="center" marginY={1}>
      <Box flexGrow={1}>
        <Typography>{hasChanges ? 'Not yet saved' : 'All changes saved'}</Typography>
      </Box>
      <Box flexGrow={1}>
        <Button
          color="primary" variant="contained"
          disabled={!hasChanges || saving} onClick={handleSave}
        >
          {saving ? 'Saving...' : 'Save'}
        </Button>
      </Box>
      <Button
        color="primary" variant="contained"
        onClick={handleConfirm}
      >
        Delete
      </Button>
    </Box>
    <Dialog
      open={openConfirmDialog}
      aria-describedby="alert-dialog-description"
    >
      <DialogContent>
        <DialogContentText id="alert-dialog-description">
          Are you sure you want to delete the item?
        </DialogContentText>
      </DialogContent>
      <DialogActions>
        <Button onClick={() => setOpenConfirmDialog(false)} autoFocus>Cancel</Button>
        <Button onClick={handleDelete}>Delete and save the changes</Button>
      </DialogActions>
    </Dialog>
  </div>
})
SectionEditor.propTypes = {
  sectionDef: PropTypes.object.isRequired,
  section: PropTypes.object,
  parent: PropTypes.any,
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
