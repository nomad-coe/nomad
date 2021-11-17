/* eslint-disable quotes */
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
import { Button } from '@material-ui/core'
import { northBase } from '../../config'

const NORTHLaunchButton = React.memo(({
  name,
  path,
  onClick,
  disabled,
  children
}) => {
  return <Button
    component="a"
    onClick={onClick}
    target="_blank"
    color="primary"
    style={{width: "6rem"}}
    variant="contained"
    disabled={disabled}
  >{children}</Button>
})

NORTHLaunchButton.propTypes = {
  name: PropTypes.string.isRequired, // The unique identitifier for this tool
  path: PropTypes.string, // Path to add to container url
  disabled: PropTypes.bool, // Whether the button is disabled
  children: PropTypes.node
}

export default NORTHLaunchButton
