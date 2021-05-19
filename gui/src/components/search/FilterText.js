
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
import { Tooltip, Typography, TextField } from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { convertSILabel } from '../../utils'
import searchQuantities from '../../searchQuantities'

const useStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    alignItems: 'center',
    margin: `${theme.spacing(1)}px 0`
  },
  rangeInput: {
    flex: '0 0 5rem',
    margin: '0 0.5rem',
    minWidth: '4rem'
  },
  dash: {
    flex: '0 0 1rem',
    textAlign: 'center'
  },
  name: {
    marginLeft: theme.spacing(1),
    minWidth: '6rem'
  }
}))
const FilterText = React.memo(({
  label,
  quantity,
  description,
  className,
  classes,
  units,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useStyles({classes: classes, theme: theme})

  // Determine the description and units
  const def = searchQuantities[quantity]
  const desc = description || def?.description || ''
  const name = label || def?.name
  let unit = def?.unit && convertSILabel(def.unit, units)

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <Tooltip title={desc}>
      <Typography className={styles.name} variant="body1">
        {`${name}:`}
      </Typography>
    </Tooltip>
    <TextField className={styles.rangeInput} margin='dense' size='small' variant='outlined' label="min"/>
    {unit && <Typography variant="body1">{unit}</Typography>}
  </div>
})

FilterText.propTypes = {
  label: PropTypes.string,
  quantity: PropTypes.string,
  description: PropTypes.string,
  className: PropTypes.string,
  classes: PropTypes.object,
  units: PropTypes.object,
  'data-testid': PropTypes.string
}

export default FilterText
