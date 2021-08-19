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
import { Tooltip, Typography } from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'

/**
 * The quantity label shown by all filter components.
 */
const useStaticStyles = makeStyles(theme => ({
  root: {
    marginBottom: theme.spacing(1)
  },
  label: {
    textTransform: 'capitalize',
    fontSize: '0.85rem',
    color: '#383838'
  }
}))
const FilterLabel = React.memo(({
  label,
  underscores,
  description,
  className,
  classes
}) => {
  const styles = useStaticStyles({classes: classes})

  // Remove underscores from name
  const finalLabel = useMemo(
    () => !underscores ? label.replaceAll('_', ' ') : label,
    [label, underscores]
  )

  return <div className={clsx(className, styles.root)}>
    <Tooltip title={description || ''} placement="bottom">
      <Typography
        className={styles.label}
        variant="button"
      >
        {finalLabel}
      </Typography>
    </Tooltip>
  </div>
})

FilterLabel.propTypes = {
  label: PropTypes.string.isRequired,
  description: PropTypes.string,
  underscores: PropTypes.bool,
  className: PropTypes.string,
  classes: PropTypes.object
}

FilterLabel.defaultProps = {
  underscores: false
}

export default FilterLabel
