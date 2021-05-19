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
import { makeStyles } from '@material-ui/core/styles'
import { Tooltip, IconButton, Button } from '@material-ui/core'
import clsx from 'clsx'

const Actions = React.memo(({
  header,
  actions,
  color,
  variant,
  size,
  justifyContent,
  className,
  classes
}) => {
  const actionsStyles = makeStyles((theme) => ({
    root: {
      display: 'flex',
      width: '100%',
      justifyContent: justifyContent,
      boxSizing: 'border-box'
    },
    spacer: {
      flexGrow: 1
    },
    iconButton: {
      marginRight: theme.spacing(1)
    }
  }))
  const styles = actionsStyles(classes)
  const buttonList = actions && actions.map((value, idx) => {
    return <Tooltip key={idx} title={value.tooltip}>
      {variant === 'icon'
        ? <IconButton
          color={color}
          size={size}
          className={styles.iconButton}
          onClick={value.onClick}
          disabled={value.disabled}
          href={value.href}
          aria-label={value.tooltip}
        >
          {value.content}
        </IconButton>
        : <Button
          color={color}
          variant={variant}
          size={size}
          className={styles.iconButton}
          onClick={value.onClick}
          disabled={value.disabled}
          href={value.href}
          aria-label={value.tooltip}
        >
          {value.content}
        </Button>
      }
    </Tooltip>
  })
  return <div className={clsx(className, styles.root)}>
    {header}
    {header && <div className={styles.spacer}></div>}
    {buttonList}
  </div>
})

Actions.propTypes = {
  header: PropTypes.any, // A text message or component to display at the left side of the actions
  actions: PropTypes.arrayOf(
    PropTypes.shape({
      content: PropTypes.any, // The content to show inside the button: text, component, icon, etc.
      href: PropTypes.string,
      tooltip: PropTypes.string,
      disabled: PropTypes.number,
      onClick: PropTypes.func
    })
  ),
  color: PropTypes.string, // The color of the MUI buttons
  variant: PropTypes.string, // The variant of the MUI buttons
  size: PropTypes.string, // Size of the MUI buttons
  justifyContent: PropTypes.string, // The flexbox justification of buttons
  className: PropTypes.string,
  classes: PropTypes.string
}

Actions.defaultProps = {
  size: 'small',
  variant: 'icon',
  justifyContent: 'flex-end'
}

export default Actions
