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
t* See the License for the specific language governing permissions and
 * limitations under the License.
 */
import React, { useState, useCallback } from 'react'
import PropTypes from 'prop-types'
import {
  Menu,
  MenuItem,
  FormControl,
  FormLabel,
  FormControlLabel,
  IconButton,
  TextField,
  Typography,
  DialogActions,
  Box,
  Tooltip,
  Radio,
  RadioGroup
} from '@material-ui/core'
import LoadingButton from './LoadingButton'
import { download } from '../../utils'
import { useErrors } from '../errors'
import { useApi } from '../api'

// TODO: The available options could be passed down to the GUI via pydantic
// model serialization.
const formats = {
  cif: {
    label: 'CIF'
  },
  xyz: {
    label: 'XYZ'
  },
  pdb: {
    label: 'PDB'
  }
}
export const wrapModes = {
  original: {
    key: 'original',
    label: 'Original',
    description: 'Original positions'
  },
  wrap: {
    key: 'wrap',
    label: 'Wrap',
    description: 'Positions are wrapped to be inside the cell respecting periodic boundary conditions'
  },
  unwrap: {
    key: 'unwrap',
    label: 'Unwrap',
    description: `Positions are reconstructed so that the structure is not split
    by periodic cell boundaries. Note that this produces meaningful results only
    if the system dimensions are smaller than the unit cell.`
  }
}

export const WrapModeRadio = ({value, onChange, disabled, className}) => {
  return <FormControl key='wrap' component="fieldset" className={className}>
    <FormLabel component="legend">Wrap mode</FormLabel>
    <RadioGroup
      value={value}
      onChange={onChange}
      >
      {Object.entries(wrapModes).map(([key, data]) =>
        <FormControlLabel
          key={key}
          value={key}
          control={<Radio color="primary" disabled={disabled}/>}
          label={<Tooltip
            title={wrapModes[key].description}>
              <span>{data.label}</span>
          </Tooltip>}
        />
      )}
    </RadioGroup>
  </FormControl>
}

WrapModeRadio.propTypes = {
  value: PropTypes.string,
  onChange: PropTypes.func,
  disabled: PropTypes.bool,
  className: PropTypes.string
}

/*
 * Menu for downloading a specific system.
 */
export const DownloadSystemMenu = React.memo(React.forwardRef(({entryId, path, anchorEl, onClose, ...rest}, ref) => {
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [format, setFormat] = useState('cif')
  const [wrapMode, setWrapMode] = useState('original')
  const [loading, setLoading] = useState(false)
  const open = Boolean(anchorEl)

  const handleChangeFormat = useCallback((event) => {
    setFormat(event.target.value)
  }, [])

  const handleChangeWrapMode = useCallback((event) => {
    setWrapMode(event.target.value)
  }, [])

  const handleClickDownload = useCallback(() => {
    setLoading(true)
    api.get(
      `systems/${entryId}?path=${path}&format=${format}&wrap_mode=${wrapMode}`,
      undefined,
      {responseType: 'blob', fullResponse: true}
    )
      .then((response) => {
        const contentDisposition = response.headers['content-disposition']
        const fileName = contentDisposition?.split('"')[1] || `system.${format}`
        download(fileName, response.data)
      })
      .finally(() => setLoading(false))
      .catch(raiseError)
  }, [entryId, path, format, wrapMode, api, raiseError])

  return <Menu
      anchorEl={anchorEl}
      getContentAnchorEl={null}
      anchorOrigin={{ vertical: 'bottom', horizontal: 'right' }}
      transformOrigin={{ vertical: 'top', horizontal: 'right' }}
      keepMounted
      open={open}
      onClose={onClose}
    >
      <Box minWidth="13rem" paddingLeft={2} paddingRight={2} paddingTop={1}>
        <Typography variant='h6' fontSize='0.9rem'>
          Download system
        </Typography>
        <Box mt={1} />
        <TextField select value={format} onChange={handleChangeFormat} size='small' variant="filled" label="Format" fullWidth>
          {Object.entries(formats).map(([key, value]) => {
            return <MenuItem key={key} value={key}>{value.label}</MenuItem>
          })}
        </TextField>
        <Box mt={2} />
        <WrapModeRadio value={wrapMode} onChange={handleChangeWrapMode} />
        <Box marginRight={-1} marginBottom={-1}>
          <DialogActions>
            <LoadingButton onClick={handleClickDownload} color="primary" loading={loading}>
              Download
            </LoadingButton>
          </DialogActions>
        </Box>
      </Box>
    </Menu>
}))

DownloadSystemMenu.propTypes = {
  entryId: PropTypes.string,
  path: PropTypes.string,
  anchorEl: PropTypes.object,
  onClose: PropTypes.func
}

/*
 * Button for downloading a specific system.
 */
export const DownloadSystemButton = React.memo(React.forwardRef(({entryId, path, ...rest}, ref) => {
  const [anchorEl, setAnchorEl] = useState(null)
  const openMenu = useCallback((event) => {
    setAnchorEl(event.currentTarget)
  }, [])
  return <>
    <IconButton size="small" onClick={openMenu} {...rest} ref={ref}/>
    <DownloadSystemMenu entryId={entryId} path={path} anchorEl={anchorEl} onClose={() => setAnchorEl(null)}/>
  </>
}))

DownloadSystemButton.propTypes = {
  entryId: PropTypes.string,
  path: PropTypes.string,
  onClose: PropTypes.func
}
