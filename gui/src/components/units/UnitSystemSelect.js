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
import React, { useCallback } from 'react'
import {FormControlLabel, RadioGroup, Radio} from '@material-ui/core'
import PropTypes from 'prop-types'
import {useUnitContext} from './UnitContext'

/**
 * Controls the unit system in the current unit context.
 */
const UnitSystemSelect = React.memo(({onChange}) => {
  const {unitSystems, selected, setSelected} = useUnitContext()

  const handleSystemChange = useCallback((event) => {
    setSelected(event.target.value)
    onChange && onChange(event)
  }, [onChange, setSelected])

  return <RadioGroup value={selected} onChange={handleSystemChange}>
    {Object.entries(unitSystems).map(([key, system]) =>
      <FormControlLabel key={key} value={key} control={<Radio />} label={system.label} />
    )}
  </RadioGroup>
})
UnitSystemSelect.propTypes = {
  onChange: PropTypes.func,
  disabled: PropTypes.bool
}

export default UnitSystemSelect
