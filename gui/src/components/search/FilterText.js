
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
import React, { useMemo } from 'react'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import { Tooltip, Typography, TextField } from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { convertSILabel } from '../../utils'
import searchQuantities from '../../searchQuantities'

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    flexDirection: 'row',
    boxSizing: 'border-box'
  },
  labelRoot: {
    fontSize: '1.1rem',
    paddingRight: '0.5rem',
    pointerEvents: 'auto',
    backgroundColor: theme.palette.background.paper,
    color: theme.palette.text.primary
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

  const labelProps = useMemo(() => ({
    shrink: true,
    classes: {
      root: styles.labelRoot
    }
  }), [styles])

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <TextField
      label={
        <Tooltip title={desc}>
          <div>{name}</div>
        </Tooltip>
      }
      variant="outlined"
      fullWidth
      InputLabelProps={labelProps}
    />
    {/* {unit && <Typography variant="body1">{unit}</Typography>} */}
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
