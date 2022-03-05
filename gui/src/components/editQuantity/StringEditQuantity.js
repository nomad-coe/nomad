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
import React, {useCallback, useState} from 'react'
import {
  Box,
  Dialog,
  DialogContent,
  DialogTitle,
  IconButton,
  InputAdornment,
  makeStyles,
  TextField
} from '@material-ui/core'
import PropTypes from 'prop-types'
import HelpOutlineIcon from '@material-ui/icons/HelpOutline'
import Markdown from '../Markdown'
import DialogActions from '@material-ui/core/DialogActions'
import Button from '@material-ui/core/Button'

const HelpDialog = React.memo(({title, description}) => {
  const [open, setOpen] = useState(false)

  return <React.Fragment>
    {description && <IconButton size="small" onClick={() => setOpen(true)}>
      {<HelpOutlineIcon fontSize='small'/>}
    </IconButton>}
    {open && <Dialog open={open}>
      <DialogTitle>
        {title}
      </DialogTitle>
      <DialogContent>
        <Markdown>{description}</Markdown>
      </DialogContent>
      <DialogActions>
        <span style={{flexGrow: 1}} />
        <Button onClick={() => setOpen(false)} color="secondary">
          Close
        </Button>
      </DialogActions>

    </Dialog>}
  </React.Fragment>
})
HelpDialog.propTypes = {
  title: PropTypes.string,
  description: PropTypes.string
}

const useHelpAdornmentStyles = makeStyles(theme => ({
  root: {},
  withOtherAdornment: {
    marginRight: theme.spacing(3)
  }
}))

const HelpAdornment = React.memo(({title, description, withOtherAdornment}) => {
  const classes = useHelpAdornmentStyles()
  return <InputAdornment
    position="end"
    className={withOtherAdornment ? classes.withOtherAdornment : classes.root}
  >
    <HelpDialog title={title} description={description}/>
  </InputAdornment>
})
HelpAdornment.propTypes = {
  withOtherAdornment: PropTypes.bool,
  title: PropTypes.string,
  description: PropTypes.string
}

const useWithHelpStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    '&:not(:hover)': {
      '& #help': {
        display: 'none'
      }
    }
  }
}))

export const WithHelp = React.memo((props) => {
  const {label, helpDescription, ...otherProps} = props
  const classes = useWithHelpStyles()
  if (!helpDescription) {
    return ''
  }
  return <Box display="flex" alignItems="center" className={classes.root}>
    <Box flexGrow={1} {...otherProps}/>
    <Box>
      <div id="help">
        <HelpDialog title={label} description={helpDescription} />
      </div>
    </Box>
  </Box>
})
WithHelp.propTypes = {
  label: PropTypes.string,
  helpDescription: PropTypes.string
}

const capitalize = (s) => {
  if (typeof s !== 'string') return ''
  return s.charAt(0).toUpperCase() + s.slice(1)
}

function getArchiveValue(quantityDef, section) {
  let value = section[quantityDef.name]
  if (value === undefined) {
    return quantityDef.default
  }
  return value
}

export function getFieldProps(quantityDef) {
  const eln = quantityDef?.m_annotations?.eln
  let name = quantityDef.name.replace(/_/g, ' ')
  let capitalizeName = capitalize(name)
  let label = (eln.length > 0 ? eln[0]?.label : undefined) || capitalizeName
  return {
    label: label,
    helpDescription: quantityDef.description
  }
}

export const TextFieldWithHelp = React.memo((props) => {
  const {withOtherAdornment, label, helpDescription, ...otherProps} = props
  const classes = useWithHelpStyles()
  return <TextField
    className={classes.root}
    InputProps={(helpDescription && {endAdornment: (
      <div id="help">
        <HelpAdornment title={label} description={helpDescription} withOtherAdornment={withOtherAdornment}/>
      </div>
    )})}
    label={label}
    {...otherProps}
  />
})
TextFieldWithHelp.propTypes = {
  withOtherAdornment: PropTypes.bool,
  label: PropTypes.string,
  helpDescription: PropTypes.string
}

export const StringEditQuantity = React.memo((props) => {
  const {quantityDef, section, onChange, ...otherProps} = props
  const value = getArchiveValue(quantityDef, section)

  const handleChange = useCallback((newValue) => {
    if (onChange) {
      onChange(newValue === '' ? undefined : newValue, section, quantityDef)
    }
  }, [onChange, quantityDef, section])

  return <TextFieldWithHelp
    fullWidth variant='filled' size='small'
    defaultValue={value !== undefined ? value : ''}
    onChange={event => handleChange(event.target.value)}
    {...getFieldProps(quantityDef)}
    {...otherProps}
  />
})
StringEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

export const StringField = React.memo((props) => {
  const {onChange, defaultValue, ...otherProps} = props

  return <TextFieldWithHelp
    fullWidth variant='filled' size='small'
    onChange={event => onChange(event.target.value)}
    defaultValue={defaultValue}
    {...otherProps}
  />
})
StringField.propTypes = {
  onChange: PropTypes.func.isRequired,
  defaultValue: PropTypes.string
}
