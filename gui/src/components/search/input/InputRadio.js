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
import { makeStyles, useTheme } from '@material-ui/core/styles'
import {
  FormControlLabel,
  RadioGroup,
  Radio,
  Tooltip,
  Typography
} from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import InputLabel from './InputLabel'
import { useFilterState } from '../FilterContext'

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    display: 'flex',
    alignItems: 'flex-start',
    justifyContent: 'center',
    flexDirection: 'column',
    boxSizing: 'border-box'
  }
}))
const InputRadio = React.memo(({
  quantity,
  label,
  description,
  options,
  className,
  classes,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useStyles({classes: classes, theme: theme})
  const [filter, setFilter] = useFilterState(quantity)

  const handleChange = useCallback((event, value) => {
    setFilter(value)
  }, [setFilter])

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <InputLabel label={label} description={description}/>
    <RadioGroup aria-label={label} name={label} value={filter || 'public'} onChange={handleChange}>
      {options && Object.entries(options).map(([key, value]) =>
        <FormControlLabel key={key} value={key} control={<Radio disabled={value.disabled}/>} label={
          <Tooltip placement="right" enterDelay={500} title={value.tooltip}>
            <Typography>{value.label}</Typography>
          </Tooltip>
        }/>
      )}
    </RadioGroup>
  </div>
})

InputRadio.propTypes = {
  quantity: PropTypes.string.isRequired,
  label: PropTypes.string,
  description: PropTypes.string,
  options: PropTypes.object, // Mapping from option name to show label and tooltip
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

export default InputRadio
