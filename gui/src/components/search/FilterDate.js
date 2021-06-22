
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
import React from 'react'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import {
  Tooltip
} from '@material-ui/core'
// import {
//   KeyboardDatePicker
// } from '@material-ui/pickers'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import searchQuantities from '../../searchQuantities'
import FilterLabel from './FilterLabel'
// import { useFilterState } from './FilterContext'

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    display: 'flex',
    alignItems: 'flex-start',
    justifyContent: 'center',
    flexDirection: 'column',
    boxSizing: 'border-box'
  },
  select: {
    width: '100%'
  },
  chips: {
    display: 'flex',
    flexWrap: 'wrap'
  },
  icon: {
    right: theme.spacing(1)
  }
}))
const FilterDate = React.memo(({
  label,
  quantity,
  description,
  visible,
  className,
  classes,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useStyles({classes: classes, theme: theme})
  // const [startDate, setStartDate] = useState(new Date())
  // const [endDate, setEndDate] = useState(new Date())
  // const [filter, setFilter] = useFilterState(quantity)
  const disabled = false

  // Determine the description and units
  const def = searchQuantities[quantity]
  const desc = description || def?.description || ''
  const title = label || def?.name

  return <Tooltip title={disabled ? 'No values available with current query.' : ''}>
    <div className={clsx(className, styles.root)} data-testid={testID}>
      <FilterLabel label={title} description={desc}/>
      {/* <KeyboardDatePicker
        autoOk
        variant="inline"
        inputVariant="outlined"
        label="With keyboard"
        format="MM/dd/yyyy"
        value={startDate}
        InputAdornmentProps={{ position: 'start' }}
        onChange={date => handleDateChange(date)}
      /> */}
    </div>
  </Tooltip>
})

FilterDate.propTypes = {
  label: PropTypes.string,
  quantity: PropTypes.string.isRequired,
  description: PropTypes.string,
  visible: PropTypes.bool.isRequired,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

export default FilterDate
