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
import { List, Divider } from '@material-ui/core'
import NorthTool from './NorthTool'
import { useTools } from './NorthPage'

const NorthLaunchButton = React.memo(function NorthLaunchButton({tools, ...props}) {
  tools = tools || []
  const northTools = useTools()
  const toolsData = Object.keys(northTools)
    .filter(tool => tools.includes(tool))
    .map(key => ({name: key, title: key, ...northTools[key]}))

  if (toolsData.length === 0) {
    return null
  }

  return (
    <List>
      {toolsData
        .map((tool, index) => (
          <div key={tool.name}>
            <NorthTool tool={tool} {...props} />
            {index !== Object.keys(toolsData).length - 1 && <Divider/>}
          </div>
        ))}
    </List>
  )
})
NorthLaunchButton.propTypes = {
  tools: PropTypes.arrayOf(PropTypes.string)
}

export default NorthLaunchButton
