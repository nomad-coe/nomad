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
import {
  Tooltip,
  IconButton,
  MenuItem,
  Select,
  Checkbox,
  FormControlLabel
} from '@material-ui/core'
import clsx from 'clsx'
import { isArray } from 'lodash'
import { useBoolState } from '../hooks'

const useActionsStyles = makeStyles((theme) => ({
  root: {
    display: 'flex',
    width: '100%',
    boxSizing: 'border-box',
    alignItems: 'center'
  }
}))
export const Actions = React.memo(({
  justifyContent,
  className,
  classes,
  children
}) => {
  const useDynamicStyles = makeStyles((theme) => ({
    root: {
      justifyContent: justifyContent
    }
  }))
  const styles = useActionsStyles({classes: classes})
  const dynamicStyles = useDynamicStyles()

  return <div className={clsx(className, styles.root, dynamicStyles.root)}>
    {children}
  </div>
})

Actions.propTypes = {
  justifyContent: PropTypes.string, // The flexbox justification of buttons
  className: PropTypes.string,
  classes: PropTypes.object,
  children: PropTypes.node
}

Actions.defaultProps = {
  justifyContent: 'flex-end'
}

const useActionHeaderStyles = makeStyles((theme) => ({
  root: {
    flexGrow: 1,
    height: '100%',
    display: 'flex',
    alignItems: 'center',
    minWidth: 0,
    marginRight: theme.spacing(0.5)
  },
  spacer: {
    flexGrow: 1,
    minWidth: 0
  }
}))
export const ActionHeader = React.memo(({
  disableSpacer,
  className,
  classes,
  children
}) => {
  const styles = useActionHeaderStyles({classes: classes})
  return <div className={clsx(className, styles.root)}>
    {children}
    {!disableSpacer && <div className={styles.spacer}></div>}
  </div>
})

ActionHeader.propTypes = {
  disableSpacer: PropTypes.bool, // Used to disable flexbox spacer
  className: PropTypes.string,
  classes: PropTypes.object,
  children: PropTypes.node
}

const useActionStyles = makeStyles((theme) => ({
  root: {
    marginRight: theme.spacing(1),
    '&:last-child': {
      marginRight: 0
    }
  }
}))
export const Action = React.memo(({
  color,
  size,
  href,
  disabled,
  onClick,
  onMouseDown,
  onMouseUp,
  tooltip,
  TooltipProps,
  className,
  classes,
  children,
  ButtonComponent,
  ButtonProps,
  'data-testid': testID
}) => {
  const styles = useActionStyles({classes: classes})

  return <Tooltip title={tooltip} {...TooltipProps}>
    <span className={clsx(className, styles.root)}>
      <ButtonComponent
        {...ButtonProps}
        color={color}
        size={size}
        onClick={onClick}
        onMouseDown={onMouseDown}
        onMouseUp={onMouseUp}
        disabled={disabled}
        href={href}
        aria-label={tooltip}
        data-testid={testID}
      >
        {children}
      </ButtonComponent>
    </span>
  </Tooltip>
})

Action.propTypes = {
  variant: PropTypes.string, // The variant of the MUI buttons
  color: PropTypes.string, // The color of the MUI buttons
  size: PropTypes.string, // Size of the MUI buttons
  href: PropTypes.string,
  disabled: PropTypes.bool,
  onClick: PropTypes.func,
  onMouseDown: PropTypes.func,
  onMouseUp: PropTypes.func,
  tooltip: PropTypes.string,
  TooltipProps: PropTypes.object,
  ButtonComponent: PropTypes.elementType,
  ButtonProps: PropTypes.object,
  className: PropTypes.string,
  classes: PropTypes.object,
  children: PropTypes.node,
  'data-testid': PropTypes.string
}

Action.defaultProps = {
  size: 'small',
  ButtonComponent: IconButton
}

/**
 * Dropdown menu action.
 */
export const ActionSelect = React.memo(({value, options, tooltip, onChange}) => {
  const [isStatsTooltipOpen, openStatsTooltip, closeStatsTooltip] = useBoolState(false)
  const items = isArray(options)
    ? Object.fromEntries(options.map(x => [x, x]))
    : options
  return <Action
    TooltipProps={{
      title: tooltip || "",
      open: isStatsTooltipOpen,
      disableHoverListener: true
  }}>
    <Select
      value={value}
      onMouseEnter={openStatsTooltip}
      onMouseLeave={closeStatsTooltip}
      onOpen={closeStatsTooltip}
      onChange={(event) => onChange && onChange(event.target.value)}
    >
      {Object.entries(items).map(([key, value]) =>
        <MenuItem key={key} value={value}>{key}</MenuItem>
      )}
    </Select>
  </Action>
})

ActionSelect.propTypes = {
  value: PropTypes.string,
  options: PropTypes.oneOfType([PropTypes.object, PropTypes.arrayOf(PropTypes.string)]),
  tooltip: PropTypes.string,
  onChange: PropTypes.func
}

/**
 * Checkbox action.
 */
export const ActionCheckbox = React.memo(({value, label, tooltip, onChange}) => {
  const styles = useActionStyles()
  return <Tooltip title={tooltip}>
      <FormControlLabel
        control={<Checkbox
          checked={value}
          onChange={(event, value) => onChange && onChange(value)}
          size="small"
        />}
        className={styles.root}
        label={label}
      />
    </Tooltip>
})

ActionCheckbox.propTypes = {
  value: PropTypes.bool,
  label: PropTypes.string,
  tooltip: PropTypes.string,
  onChange: PropTypes.func
}
