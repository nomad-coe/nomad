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
import PropTypes from 'prop-types'
import { Tooltip, IconButton, Button, Box, makeStyles } from '@material-ui/core'
import clsx from 'clsx'

export default function Actions({actions, color, variant, size, justifyContent, className, classes}) {
  const actionsStyles = makeStyles((theme) => ({
    root: {
      display: 'flex',
      width: '100%',
      justifyContent: justifyContent
    },
    iconButton: {
      marginRight: theme.spacing(1)
    }
  }))
  const styles = actionsStyles(classes)
  const buttonList = actions.map((value, idx) => {
    return <Tooltip key={idx} title={value.tooltip}>
      {variant === 'icon'
        ? <IconButton
          color={color}
          size={size}
          className={styles.iconButton}
          onClick={value.onClick}
          disabled={value.disabled}
          href={value.href}>
          {value.content}
        </IconButton>
        : <Button
          color={color}
          variant={variant}
          size={size}
          className={styles.iconButton}
          onClick={value.onClick}
          disabled={value.disabled}
          href={value.href}>
          {value.content}
        </Button>
      }
    </Tooltip>
  })
  return <Box className={clsx(className, styles.root)}>
    {buttonList}
  </Box>
}

Actions.propTypes = {
  actions: PropTypes.array,
  color: PropTypes.string,
  variant: PropTypes.string,
  size: PropTypes.string,
  justifyContent: PropTypes.string,
  className: PropTypes.string,
  classes: PropTypes.string
}

Actions.defaultProps = {
  size: 'small',
  variant: 'icon',
  justifyContent: 'flex-end'
}
