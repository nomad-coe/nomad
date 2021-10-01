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
import { Button, Tooltip } from '@material-ui/core'

/** A button wrapped in a tooltip. All further props are forwarded to the button component. */
const TooltipButton = React.forwardRef(function TooltipButton(props, ref) {
  const {title, component, ...moreProps} = props
  const button = React.createElement(component || Button, {...moreProps, ref: ref})
  return <Tooltip title={title}>
    {button}
  </Tooltip>
})
TooltipButton.propTypes = {
  /** The tooltip title. */
  title: PropTypes.string.isRequired,
  /** The button component. Default is MUI Button. */
  component: PropTypes.elementType,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
}

export default TooltipButton
