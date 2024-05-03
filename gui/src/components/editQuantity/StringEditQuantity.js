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
import LaunchIcon from '@material-ui/icons/Launch'
import {getDisplayLabel} from "../../utils"
import {useRecoilValue} from "recoil"
import {configState} from "../archive/ArchiveBrowser"

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

export const HelpAdornment = React.memo(({title, description, withOtherAdornment}) => {
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
  return <Box display="flex" alignItems="center" className={classes.root}>
    <Box flexGrow={1} {...otherProps}/>
    {helpDescription && (
      <Box>
        <div id="help">
          <HelpDialog title={label} description={helpDescription} />
        </div>
      </Box>
    )}
  </Box>
})
WithHelp.propTypes = {
  label: PropTypes.string,
  helpDescription: PropTypes.string
}

export function getFieldProps(quantityDef) {
  const eln = quantityDef?.m_annotations?.eln
  const {component, ...elnProps} = eln?.[0] || {}
  elnProps.unit = elnProps.unit || quantityDef.unit
  return {
    helpDescription: quantityDef.description,
    ...elnProps
  }
}

export const TextFieldWithHelp = React.memo(React.forwardRef((props, ref) => {
  const {withOtherAdornment, label, helpDescription, 'data-testid': TestId, ...otherProps} = props
  const classes = useWithHelpStyles()
  return <TextField
    inputRef={ref}
    className={classes.root}
    InputProps={(helpDescription && {endAdornment: (
      <div id="help">
        <HelpAdornment title={label} description={helpDescription} withOtherAdornment={withOtherAdornment}/>
      </div>
    )})}
    label={label}
    data-testid={TestId}
    {...otherProps}
    value={otherProps.value || ''}
  />
}))
TextFieldWithHelp.propTypes = {
  withOtherAdornment: PropTypes.bool,
  label: PropTypes.string,
  value: PropTypes.string,
  helpDescription: PropTypes.string,
  'data-testid': PropTypes.string
}

export const StringEditQuantity = React.memo((props) => {
  const {quantityDef, onChange, ...otherProps} = props
  const config = useRecoilValue(configState)
  const label = getDisplayLabel(quantityDef, true, config?.showMeta)

  const handleChange = useCallback((value) => {
    if (onChange) {
      onChange(value === '' ? undefined : value)
    }
  }, [onChange])
  return <TextFieldWithHelp
    fullWidth variant='filled' size='small'
    onChange={event => handleChange(event.target.value)}
    {...getFieldProps(quantityDef)}
    {...otherProps}
    label={label}
  />
})
StringEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.string,
  onChange: PropTypes.func
}

export const StringField = React.memo((props) => {
  const {onChange, ...otherProps} = props

  const handleChange = useCallback(event => {
    const value = event.target.value
    if (onChange) {
      onChange(value === '' ? undefined : value)
    }
  }, [onChange])

  return <TextFieldWithHelp
    fullWidth variant='filled' size='small'
    onChange={handleChange}
    {...otherProps}
  />
})
StringField.propTypes = {
  onChange: PropTypes.func
}

export const TextFieldWithLinkButton = React.memo(React.forwardRef((props, ref) => {
  const {withOtherAdornment, label, value, helpDescription, 'data-testid': TestId, ...otherProps} = props
  const classes = useWithHelpStyles()

  const validateURL = useCallback((value) => {
    try {
      return Boolean(new URL(value))
    } catch (e) {
      return false
    }
  }, [])

  return <TextField
    defaultValue={value}
    error={value !== undefined && !validateURL(value)}
    helperText={value === undefined || validateURL(value) ? '' : 'Invalid URL string!'}
    inputRef={ref}
    className={classes.root}
    InputProps={{endAdornment:
      <>
        { helpDescription &&
        <div id="help">
          <HelpAdornment title={label} description={helpDescription} withOtherAdornment={withOtherAdornment}/>
        </div>
        }
        {
          validateURL(value) &&
          <IconButton aria-label="open_link" onClick={() => window.open(value, '_blank')}>
            <LaunchIcon />
          </IconButton>
        }
      </>
    }}
    label={label}
    data-testid={'URLEditQuantity'}
    {...otherProps}
  />
}))
TextFieldWithLinkButton.propTypes = {
  withOtherAdornment: PropTypes.bool,
  label: PropTypes.string,
  value: PropTypes.string,
  helpDescription: PropTypes.string,
  'data-testid': PropTypes.string
}

export const URLEditQuantity = React.memo((props) => {
  const {quantityDef, onChange, ...otherProps} = props
  const config = useRecoilValue(configState)
  const label = getDisplayLabel(quantityDef, true, config?.showMeta)

  const handleChange = useCallback((value) => {
    if (onChange) {
      onChange(value === '' ? undefined : value)
    }
  }, [onChange])

  return <TextFieldWithLinkButton
    fullWidth variant='filled' size='small'
    onChange={event => handleChange(event.target.value)}
    {...getFieldProps(quantityDef)}
    {...otherProps}
    label={label}
  />
})
URLEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.string,
  onChange: PropTypes.func
}
