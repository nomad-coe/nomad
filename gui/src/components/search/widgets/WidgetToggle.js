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
import { cloneDeep } from 'lodash'
import AddCircleIcon from '@material-ui/icons/AddCircle'
import CancelIcon from '@material-ui/icons/Cancel'
import { useSearchContext } from '../SearchContext'
import { Action } from '../../Actions'

/**
 * Toggle for hiding/showing a particular widget.
 */
const WidgetToggle = React.memo(({quantity, disabled, 'data-testid': testID}) => {
  const { useWidgetValue, useAddWidget, useRemoveWidget, filterData } = useSearchContext()
  const widget = useWidgetValue(quantity)
  const addWidget = useAddWidget()
  const removeWidget = useRemoveWidget()
  const hasWidget = !!widget
  const widgetDefault = filterData[quantity]?.widget

  // If there is no default widget configured, then do not show the toggle.
  if (!widgetDefault) return null

  return <Action
    tooltip={hasWidget
      ? 'Remove widget from the dashboard.'
      : 'Add a default widget for this quantity in the dashboard.'
    }
    disabled={disabled}
    data-testid={testID}
    onClick={() => {
      const config = {
        id: quantity,
        editing: false,
        visible: true,
        quantity: quantity,
        ...cloneDeep(widgetDefault),
        scale: filterData[quantity].scale
      }
      if (hasWidget) {
        removeWidget(quantity)
      } else {
        addWidget(quantity, config)
      }
    }}
  >
    {hasWidget
      ? <CancelIcon fontSize="small"/>
      : <AddCircleIcon fontSize="small"/>}
  </Action>
})

WidgetToggle.propTypes = {
  quantity: PropTypes.string.isRequired,
  disabled: PropTypes.bool,
  'data-testid': PropTypes.string
}

export default WidgetToggle
