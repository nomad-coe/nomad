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
import React, { useMemo, useState } from 'react'
import { makeStyles } from '@material-ui/core/styles'
import { Tooltip, Typography, Select, MenuItem } from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { Actions } from '../../Actions'
import { useSearchContext } from '../SearchContext'

/**
 * The quantity label shown by all filter components.
 */
const useStaticStyles = makeStyles(theme => ({
  root: {
    marginBottom: theme.spacing(0.5),
    height: '2.5rem',
    width: '100%'
  },
  label: {
    textTransform: 'capitalize',
    fontSize: '0.9rem',
    color: '#383838'
  }
}))
const scales = {
  'linear': 1,
  '1/2': 0.5,
  '1/4': 0.25,
  '1/8': 0.125
}
const FilterLabel = React.memo(({
  label,
  underscores,
  description,
  disableScale,
  scale,
  onChangeScale,
  className,
  classes
}) => {
  const styles = useStaticStyles({classes: classes})
  const {useIsStatisticsEnabled} = useSearchContext()
  const isStatisticsEnabled = useIsStatisticsEnabled()
  const [open, setOpen] = useState(false)

  // Remove underscores from name
  const finalLabel = useMemo(
    () => !underscores ? label.replace(/_/g, ' ') : label,
    [label, underscores]
  )

  // The tooltip needs to be controlled: otherwise it won't close as we open the
  // select menu

  return <Actions
    className={clsx(className, styles.root)}
    header={
      <Tooltip title={description || ''} placement="bottom">
        <Typography
          className={styles.label}
          variant="button"
        >
          {finalLabel}
        </Typography>
      </Tooltip>
    }
  >
    {(!disableScale && isStatisticsEnabled) &&
      <Tooltip open={open} title="Select the scaling of the statistics.">
        <Select
          value={scale}
          onChange={(event) => onChangeScale(event.target.value)}
          onMouseEnter={() => setOpen(true)}
          onMouseLeave={() => setOpen(false)}
          onOpen={() => setOpen(false)}
          displayEmpty
          name="scale power"
        >
          {Object.entries(scales).map(([key, value]) => (
            <MenuItem key={key} value={value}>{key}</MenuItem>
          ))}
        </Select>
      </Tooltip>
    }
  </Actions>
})

FilterLabel.propTypes = {
  label: PropTypes.string.isRequired,
  description: PropTypes.string,
  underscores: PropTypes.bool,
  disableScale: PropTypes.bool,
  scale: PropTypes.oneOf(Object.values(scales)),
  onChangeScale: PropTypes.func,
  className: PropTypes.string,
  classes: PropTypes.object
}

FilterLabel.defaultProps = {
  underscores: false,
  disableScale: false,
  scale: 1
}

export default FilterLabel
