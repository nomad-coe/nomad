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
import { Chip } from '@material-ui/core'
import PropTypes from 'prop-types'

/**
 * Thin wrapper for MUI Chip that is used for displaying (and possibly removing)
 * filter values.
 */
const useStyles = makeStyles(theme => ({
  root: {
    padding: theme.spacing(0.5)
  }
}))
const FilterChip = React.memo(({
  label,
  onDelete,
  className,
  locked
}) => {
  const styles = useStyles()

  return <div className={clsx(className, styles.root)}>
    <Chip
      label={label}
      onDelete={locked ? undefined : onDelete}
      color={'primary'}
      className={styles.root}
      icon={locked ? <LockIcon/> : undefined}
    />
  </div>
})

FilterChip.propTypes = {
  label: PropTypes.oneOfType([PropTypes.string, PropTypes.number]),
  onDelete: PropTypes.func,
  className: PropTypes.string,
  locked: PropTypes.bool
}

export default FilterChip
