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
import React, { useCallback, useEffect, useState } from 'react'
import PropTypes from 'prop-types'
import { makeStyles, Typography, Button, TextField } from '@material-ui/core'
import EditIcon from '@material-ui/icons/Edit'
import WithButton from '../utils/WithButton'

/**
 * Displays the upload name and controls for updating it.
 */
const useStyles = makeStyles(theme => ({
  row: {
    display: 'flex',
    alignItems: 'center',
    '& :not(:first-child)': {
      marginLeft: theme.spacing(2)
    }
  }
}))

export const UploadName = React.memo(({upload_name, onChange}) => {
  const [edit, setEdit] = useState(false)
  const [value, setValue] = useState()
  const styles = useStyles()

  // Any changes to the upload name through the prop override the shown value.
  useEffect(() => {
    setValue(upload_name)
  }, [upload_name])

  // The upload name is trimmed, and an untouched input field is considered as
  // an empty string. Currently the only way to unset upload name is to post an
  // empty string, posting 'undefined' gives an error.
  const handleSave = useCallback(() => {
    setEdit(false)
    const trimmedValue = value?.trim() || ''
    onChange && onChange(trimmedValue)
    setValue(trimmedValue)
  }, [onChange, value])

  const handleEdit = useCallback(() => {
    setEdit(true)
  }, [])

  const handleChange = useCallback((event) => {
    setValue(event.target.value)
  }, [])

  return edit
    ? <div className={styles.row}>
        <TextField value={value || ""} onChange={handleChange} fullWidth/>
        <Button size="small" variant="contained" onClick={handleSave}>save</Button>
      </div>
    : <WithButton size="small" icon={<EditIcon/>} onClick={handleEdit}>
        <Typography variant="h6">
          {value || <i>unnamed upload</i>}
        </Typography>
      </WithButton>
})
UploadName.propTypes = {
  upload_name: PropTypes.string,
  onChange: PropTypes.func
}

export default UploadName
