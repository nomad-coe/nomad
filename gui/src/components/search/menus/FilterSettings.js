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
import { makeStyles } from '@material-ui/core/styles'
import {
  Checkbox,
  MenuItem,
  FormControlLabel
} from '@material-ui/core'
import { ToggleButton, ToggleButtonGroup } from '@material-ui/lab'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { useSearchContext } from '../SearchContext'

/**
 * Menu showing the Filter settings.
 */
const useStyles = makeStyles((theme) => {
  return {
    root: {
    },
    menuItem: {
      width: '10rem'
    },
    systems: {
      margin: theme.spacing(2),
      marginTop: theme.spacing(1)
    }
  }
})
const FilterSettings = React.memo(({
  className,
  classes,
  'data-testid': testID
}) => {
  const styles = useStyles(classes)
  const {
    useIsStatisticsEnabled,
    useSetIsStatisticsEnabled,
    useStatisticsCountMode,
    useSetStatisticsCountMode
  } = useSearchContext()
  const [isStatisticsEnabled, setIsStatisticsEnabled] = [useIsStatisticsEnabled(), useSetIsStatisticsEnabled()]
  const [statisticsCountMode, setStatisticsCountMode] = [useStatisticsCountMode(), useSetStatisticsCountMode()]

  const handleStatsChange = useCallback((event, value) => {
    setIsStatisticsEnabled(value)
  }, [setIsStatisticsEnabled])

  return <div className={clsx(styles.root, className)} data-testid={testID}>
    <MenuItem>
      <FormControlLabel
        control={<Checkbox
          checked={isStatisticsEnabled}
          onChange={handleStatsChange}
        />}
        label="Show statistics"
      />
    </MenuItem>
    <MenuItem>
      <ToggleButtonGroup
        size="small"
        exclusive
        value={statisticsCountMode}
        onChange={(event, value) => { if (value) setStatisticsCountMode(value) }}
      >
        <ToggleButton value="tooltip">Tooltip</ToggleButton>
        <ToggleButton value="fixed">Fixed</ToggleButton>
        <ToggleButton value="none">None</ToggleButton>
      </ToggleButtonGroup>
    </MenuItem>
  </div>
})

FilterSettings.propTypes = {
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

export default FilterSettings
