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
import JupyterIcon from '@material-ui/icons/Settings'
import { IconButton, Tooltip } from '@material-ui/core'
import { northBase } from '../../config'

const NorthButton = React.forwardRef(function NorthButton(props, ref) {
  const {path, anonymous, component, ...moreProps} = props
  if (!northBase) {
    return ''
  }

  if (!moreProps.children) {
    if ((moreProps.component || IconButton) === IconButton) {
      moreProps.children = <Tooltip title="Open with Jupyter">
        <JupyterIcon />
      </Tooltip>
    } else {
      moreProps.children = 'Open with Jupyter'
    }
  }

  const redirectUrl = `${northBase}/hub/user-redirect/${path || ''}`
  const url = anonymous
    ? `${northBase}/hub/anonymous_login?next=${encodeURIComponent(redirectUrl)}`
    : redirectUrl

  return React.createElement(component || IconButton, {
    target: '_blank',
    component: 'a',
    href: url,
    ...moreProps,
    ref: ref
  })
})
NorthButton.propTypes = {
  path: PropTypes.string,
  anonymous: PropTypes.bool,
  component: PropTypes.elementType,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
}

export default NorthButton
