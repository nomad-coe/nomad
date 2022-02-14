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
import React, {useCallback, useEffect, useState} from 'react'
import {TextField, makeStyles, InputAdornment} from '@material-ui/core'
import PropTypes from 'prop-types'

const useStyles = makeStyles(theme => ({
  editQantity: {
    display: 'block',
    width: '100%'
  },
  adornment: {
    marginRight: theme.spacing(3)
  }
}))

export const StringEditQantity = React.memo((props) => {
  const classes = useStyles()
  const {quantityDef, section, multiline, onChange} = props
  const [value, setValue] = useState()

  useEffect(() => {
    setValue(section[quantityDef.key])
  }, [quantityDef.key, section])

  const handleChange = useCallback((value) => {
    setValue(value)
    if (onChange) {
      onChange(value, section, quantityDef)
    }
  }, [onChange, quantityDef, section])

  return <TextField
    fullWidth variant="filled" size='small' value={value || ''} label={quantityDef?.name} multiline={multiline} minRows={4}
    InputProps={{endAdornment: <InputAdornment className={classes.adornment} position='end'>{quantityDef?.unit}</InputAdornment>}}
    placeholder={quantityDef?.description}
    onChange={event => handleChange(event.target.value)} style={{pointerEvents: 'auto'}}>
  </TextField>
})
StringEditQantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func,
  multiline: PropTypes.bool
}
