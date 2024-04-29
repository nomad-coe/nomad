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

import React, { useCallback, useEffect, useRef, useState } from 'react'
import PropTypes from 'prop-types'
import { TextField } from '@material-ui/core'
import { isNil, isEmpty } from 'lodash'
import YAML from 'yaml'

/**
 * Form used for editing a YAML or JSON config.
 */
const InputConfig = React.memo(({data, format, maxRows, minRows, onChange, error, onError, readOnly}) => {
  const controlledError = useRef(error !== undefined)
  const [serialized, setSerialized] = useState()
  const [errorInternal, setErrorInternal] = useState()
  const errorFinal = error || errorInternal
  const formatLib = {
    'JSON': JSON,
    'YAML': YAML
  }[format]

  useEffect(() => {
    setSerialized(isEmpty(data)
      ? ''
      : formatLib.stringify(
        data,
        (k, v) => { return v === Infinity ? "Infinity" : v },
        2)
    )
  }, [data, formatLib])

  const handleError = useCallback((e) => {
    if (!controlledError.current) {
      setErrorInternal(e)
      onError && onError(e)
    } else {
      if (!isNil(e)) onError && onError(e)
    }
  }, [onError])

  const handleChange = useCallback((event) => {
    const value = event.target.value
    setSerialized(value)
    try {
      const data = formatLib.parse(value)
      onChange && onChange(data)
      handleError(null)
    } catch (e) {
      handleError(`This is not ${format}: ` + e)
    }
  }, [onChange, handleError, formatLib, format])

  return (
    <TextField
      fullWidth
      label={format}
      error={!!errorFinal}
      helperText={errorFinal}
      variant="filled"
      multiline
      maxRows={maxRows}
      minRows={minRows}
      value={serialized}
      onChange={handleChange}
      InputProps={{
        readOnly: readOnly
      }}
    />
  )
})
InputConfig.propTypes = {
  data: PropTypes.oneOfType([PropTypes.object, PropTypes.array]),
  format: PropTypes.oneOf(['JSON', 'YAML']),
  maxRows: PropTypes.number,
  minRows: PropTypes.number,
  onChange: PropTypes.func,
  error: PropTypes.string,
  onError: PropTypes.func,
  readOnly: PropTypes.bool
}

InputConfig.defaultProps = {
  format: 'JSON',
  maxRows: 20
}

export default InputConfig
