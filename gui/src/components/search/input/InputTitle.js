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
import { useSearchContext } from '../SearchContext'
import { inputSectionContext } from './InputSection'
import { useUnits, Unit } from '../../../units'

/**
 * Title for a search filter. The stylized name and description are
 * automatically retrieved from metainfo or the filter data.
 */
const useStaticStyles = makeStyles(theme => ({
  root: {
  },
  title: {
    fontWeight: 600,
    color: theme.palette.grey[800]
  }
}))
const InputTitle = React.memo(({
  quantity,
  description,
  variant,
  TooltipProps,
  onMouseDown,
  onMouseUp,
  anchored,
  className,
  classes,
  style
}) => {
  const styles = useStaticStyles({classes: classes})
  const { filterData } = useSearchContext()
  const sectionContext = useContext(inputSectionContext)
  const units = useUnits()
  const section = sectionContext?.section

  // Remove underscores from name
  const finalLabel = useMemo(() => {
    const prefix = section && anchored ? `${filterData[section]?.label} ` : ''
    let label = `${prefix}${filterData[quantity]?.label}`
    const unit = filterData[quantity]?.unit
    if (unit) {
      const unitDef = new Unit(unit)
      label = `${label} (${unitDef.toSystem(units).label()})`
    }
    return label
  }, [filterData, quantity, units, section, anchored])

  const finalDescription = description || filterData[quantity].description

  return <Tooltip title={finalDescription || ''} placement="bottom" {...(TooltipProps || {})}>
    <Typography
      noWrap
      className={clsx(className, styles.root, (!section || anchored) && styles.title)}
      variant={variant}
      onMouseDown={onMouseDown}
      onMouseUp={onMouseUp}
      style={style}
    >
      {finalLabel}
    </Typography>
  </Tooltip>
})

InputTitle.propTypes = {
  quantity: PropTypes.string.isRequired,
  description: PropTypes.string,
  variant: PropTypes.string,
  anchored: PropTypes.bool,
  className: PropTypes.string,
  classes: PropTypes.object,
  style: PropTypes.object,
  TooltipProps: PropTypes.object, // Properties forwarded to the Tooltip
  onMouseDown: PropTypes.func,
  onMouseUp: PropTypes.func
}

InputTitle.defaultProps = {
  variant: 'body2'
}

export default InputTitle
