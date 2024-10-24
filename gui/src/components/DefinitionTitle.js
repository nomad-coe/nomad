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
import { makeStyles } from '@material-ui/core/styles'
import { Typography, Tooltip } from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import Ellipsis from './visualization/Ellipsis'

/**
 * Simple component for displaying titles and corresponding tooltips for
 * metainfo definitions.
 */
const useStaticStyles = makeStyles(theme => ({
  root: {
  },
  text: {
  },
  subtitle2: {
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
export const DefinitionTitle = React.memo(({
  label,
  description,
  variant,
  TooltipProps,
  onMouseDown,
  onMouseUp,
  className,
  classes,
  rotation,
  section,
  noWrap
}) => {
  const styles = useStaticStyles({classes})

  return <Tooltip title={description || ''} interactive enterDelay={400} enterNextDelay={400} {...(TooltipProps || {})}>
    <div className={clsx(className, styles.root,
      rotation === 'right' && styles.right,
      rotation === 'down' && styles.down,
      rotation === 'up' && styles.up
    )}>
      <Typography
        noWrap={noWrap}
        className={clsx(styles.text, (!section) && (variant === "subtitle2") && styles.subtitle2)}
        variant={variant}
        onMouseDown={onMouseDown}
        onMouseUp={onMouseUp}
      >
        <Ellipsis>{label}</Ellipsis>
      </Typography>
    </div>
  </Tooltip>
})

DefinitionTitle.propTypes = {
  quantity: PropTypes.string,
  label: PropTypes.string,
  description: PropTypes.oneOfType([PropTypes.string, PropTypes.node]),
  variant: PropTypes.string,
  className: PropTypes.string,
  classes: PropTypes.object,
  rotation: PropTypes.oneOf(['up', 'right', 'down']),
  TooltipProps: PropTypes.object, // Properties forwarded to the Tooltip
  onMouseDown: PropTypes.func,
  onMouseUp: PropTypes.func,
  placement: PropTypes.string,
  noWrap: PropTypes.bool,
  section: PropTypes.string
}

DefinitionTitle.defaultProps = {
  variant: 'body2',
  rotation: 'right',
  noWrap: true
}
