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
import React, { useMemo, useContext } from 'react'
import { makeStyles } from '@material-ui/core/styles'
import { Typography, Tooltip } from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { useSearchContext } from './SearchContext'
import { inputSectionContext } from './input/InputSection'
import { Unit } from '../units/Unit'
import { useUnitContext } from '../units/UnitContext'

/**
 * Title for a metainfo quantity or section that is used in a search context.
 * By default the label, description and unit are automatically retrieved from
 * the filter config.
 */
const useStaticStyles = makeStyles(theme => ({
  root: {
  },
  title: {
    fontWeight: 600,
    color: theme.palette.grey[800]
  },
  right: {
    overflow: 'hidden'
  },
  down: {
    overflow: 'hidden',
    writingMode: 'vertical-rl',
    textOrientation: 'mixed'
  },
  up: {
    overflow: 'hidden',
    writingMode: 'vertical-rl',
    textOrientation: 'mixed',
    transform: 'rotate(-180deg)'
  }
}))
const FilterTitle = React.memo(({
  quantity,
  label,
  description,
  unit,
  variant,
  TooltipProps,
  onMouseDown,
  onMouseUp,
  className,
  classes,
  rotation,
  disableUnit
}) => {
  const styles = useStaticStyles({classes: classes})
  const { filterData } = useSearchContext()
  const sectionContext = useContext(inputSectionContext)
  const {units} = useUnitContext()
  const section = sectionContext?.section

  // Create the final label
  const finalLabel = useMemo(() => {
    let finalLabel = label || filterData[quantity]?.label
    if (!disableUnit) {
      let finalUnit
      if (unit) {
        finalUnit = new Unit(unit).label()
      } else if (filterData[quantity]?.unit) {
        finalUnit = new Unit(filterData[quantity].unit).toSystem(units).label()
      }
      if (finalUnit) {
        finalLabel = `${finalLabel} (${finalUnit})`
      }
    }
    return finalLabel
  }, [filterData, quantity, units, label, unit, disableUnit])

  // Determine the final description
  const finalDescription = description || filterData[quantity]?.description || ''

  return <Tooltip title={finalDescription} placement="bottom" {...(TooltipProps || {})}>
    <div className={clsx(
      rotation === 'right' && styles.right,
      rotation === 'down' && styles.down,
      rotation === 'up' && styles.up
    )}>
      <Typography
        noWrap
        className={clsx(className, styles.root, (!section) && styles.title)}
        variant={variant}
        onMouseDown={onMouseDown}
        onMouseUp={onMouseUp}
      >
        {finalLabel}
      </Typography>
    </div>
  </Tooltip>
})

FilterTitle.propTypes = {
  quantity: PropTypes.string,
  label: PropTypes.string,
  unit: PropTypes.string,
  description: PropTypes.string,
  variant: PropTypes.string,
  className: PropTypes.string,
  classes: PropTypes.object,
  rotation: PropTypes.oneOf(['up', 'right', 'down']),
  disableUnit: PropTypes.bool,
  TooltipProps: PropTypes.object, // Properties forwarded to the Tooltip
  onMouseDown: PropTypes.func,
  onMouseUp: PropTypes.func
}

FilterTitle.defaultProps = {
  variant: 'body2',
  rotation: 'right'
}

export default FilterTitle
