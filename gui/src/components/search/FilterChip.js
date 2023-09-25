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
import clsx from 'clsx'
import { makeStyles } from '@material-ui/core/styles'
import LockIcon from '@material-ui/icons/Lock'
import { Chip, Tooltip, Typography } from '@material-ui/core'
import PropTypes from 'prop-types'
import FilterTitle from './FilterTitle'

/**
 * Thin wrapper for MUI Chip that is used for displaying (and possibly removing)
 * filter values.
 */
const useStyles = makeStyles(theme => ({
  root: {
    padding: theme.spacing(0.5),
    boxSizing: 'border-box',
    maxWidth: '100%'
  },
  chip: {
    padding: theme.spacing(0.5),
    maxWidth: '100%'
  }
}))
export const FilterChip = React.memo(({
  label,
  onDelete,
  color,
  className,
  locked
}) => {
  const styles = useStyles()

  return <div className={clsx(className, styles.root)}>
    <Tooltip title={locked ? 'This filter is locked in the current search context: it cannot be removed or modified.' : ''}>
      <Chip
        label={label}
        onDelete={locked ? undefined : onDelete}
        color={locked ? undefined : color}
        className={styles.chip}
        icon={locked ? <LockIcon/> : undefined}
      />
    </Tooltip>
  </div>
})

FilterChip.propTypes = {
  label: PropTypes.oneOfType([PropTypes.string, PropTypes.number]),
  onDelete: PropTypes.func,
  color: PropTypes.string,
  className: PropTypes.string,
  locked: PropTypes.bool
}
FilterChip.defaultProps = {
  color: 'primary'
}

export default FilterChip

/**
 * Used to group several related filter chips inside one container.
 */
const useFilterChipGroupStyles = makeStyles(theme => ({
  root: {
    borderRadius: theme.spacing(1),
    boxSizing: 'border-box',
    width: '100%'
  },
  paper: {
    display: 'flex',
    flexDirection: 'row',
    flexWrap: 'wrap',
    alignItems: 'center',
    width: '100%',
    boxSizing: 'border-box'
  },
  title: {
    marginLeft: theme.spacing(0.5),
    color: theme.palette.grey[600]
  }
}))
export const FilterChipGroup = React.memo(({
  quantity,
  className,
  children
}) => {
  const styles = useFilterChipGroupStyles()

  return <div className={clsx(className, styles.root)}>
    <FilterTitle
      quantity={quantity}
      variant="caption"
      className={styles.title}
    />
    <div className={styles.paper}>
      {children}
    </div>
  </div>
})

FilterChipGroup.propTypes = {
  quantity: PropTypes.string,
  color: PropTypes.string,
  className: PropTypes.string,
  children: PropTypes.node
}
FilterChipGroup.defaultProps = {
  color: 'primary'
}

/**
 * Operators between filter chips.
 */
const useFilterOpStyles = makeStyles(theme => ({
  root: {
    fontSize: '0.65rem'
  }
}))
const FilterOp = React.memo(({className, children}) => {
  const styles = useFilterOpStyles()
  return <Typography variant="caption" className={clsx(className, styles.root)}>{children}</Typography>
})

FilterOp.propTypes = {
  className: PropTypes.string,
  children: PropTypes.node
}

export const FilterAnd = React.memo(() => {
  return <FilterOp>AND</FilterOp>
})

export const FilterOr = React.memo(() => {
  return <FilterOp>OR</FilterOp>
})
