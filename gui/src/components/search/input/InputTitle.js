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
import { makeStyles } from '@material-ui/core/styles'
import { Typography, Tooltip } from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import searchQuantities from '../../../searchQuantities'
import { useSearchContext } from '../SearchContext'

/**
 * Title for a search filter. The stylized name and description are
 * automatically retrieved from metainfo or the filter data.
 */
const useStaticStyles = makeStyles(theme => ({
  root: {
  },
  capitalize: {
    textTransform: 'capitalize'
  }
}))
const InputTitle = React.memo(({
  quantity,
  description,
  variant,
  underscores,
  capitalize,
  TooltipProps,
  onMouseDown,
  onMouseUp,
  className,
  classes,
  style
}) => {
  const styles = useStaticStyles({classes: classes})
  const { filterData } = useSearchContext()

  // Remove underscores from name
  const finalLabel = useMemo(() => {
    let label = filterData[quantity]?.label
    if (!label) {
      label = searchQuantities[quantity]?.name || quantity
      label = !underscores ? label.replace(/_/g, ' ') : label
    }
    return label
  }, [filterData, quantity, underscores])

  const finalDescription = description || searchQuantities[quantity]?.description

  return <Tooltip title={finalDescription || ''} placement="bottom" {...(TooltipProps || {})}>
    <Typography
      noWrap
      className={clsx(className, styles.root, capitalize && styles.capitalize)}
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
  underscores: PropTypes.bool,
  className: PropTypes.string,
  classes: PropTypes.object,
  style: PropTypes.object,
  capitalize: PropTypes.bool,
  TooltipProps: PropTypes.object, // Properties forwarded to the Tooltip
  onMouseDown: PropTypes.func,
  onMouseUp: PropTypes.func
}

InputTitle.defaultProps = {
  capitalize: true,
  variant: 'body2'
}

export default InputTitle
