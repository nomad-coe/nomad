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
import React, { useEffect, useState } from 'react'
import PropTypes from 'prop-types'
import { Box } from '@material-ui/core'
import ElectronicStructureOverview from './ElectronicStructureOverview'
import { useRouteMatch, Route } from 'react-router-dom'
import { withApi } from '../api'

const MaterialPageContent = withApi(false, true)(({fixed, api, materialId, raiseError}) => {
  const props = fixed ? {maxWidth: 1200} : {}
  const [data, setData] = useState()

  // Load the data parallelly from API on first render
  useEffect(() => {
    Promise.all([
      api.encyclopediaBasic(materialId),
      api.encyclopediaCalculations(materialId)
    ]).then((results) => {
      setData({
        basic: results[0],
        calculations: results[1]
      })
    }).catch(error => {
      if (error.name === 'DoesNotExist') {
        raiseError(error)
      }
    })
  }, [api, materialId, raiseError])

  return <Box padding={3} margin="auto" {...props}>
    <ElectronicStructureOverview data={data}></ElectronicStructureOverview>
  </Box>
})
MaterialPageContent.propTypes = ({
  materialId: PropTypes.string,
  api: PropTypes.func,
  raiseError: PropTypes.func,
  fixed: PropTypes.bool
})
function MaterialPage() {
  const { path } = useRouteMatch()

  return (
    <Route
      path={`${path}/:materialId?/:tab?`}
      render={({match: {params: {materialId, tab = 'overview'}}}) => {
        if (materialId) {
          return (
            <React.Fragment>
              <MaterialPageContent fixed={true} materialId={materialId}>
              </MaterialPageContent>
            </React.Fragment>
          )
        } else {
          return ''
        }
      }}
    />
  )
}

export default MaterialPage
