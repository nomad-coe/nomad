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
import Grid from '@material-ui/core/Grid'
import { Quantity } from '../QuantityHistogram'
import { withStyles } from '@material-ui/core'
import { searchContext } from '../SearchContext'
import { compose } from 'recompose'
import { withApi } from '../../api'

class UploadersList extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired
  }

  static styles = theme => ({
    root: {
      marginTop: theme.spacing(2)
    }
  })

  static contextType = searchContext

  render() {
    const {state: {usedMetric}} = this.context

    return (
      <Grid>
        <Quantity quantity="origin" title="Uploaders/origin" scale={1} metric={usedMetric} />
      </Grid>
    )
  }
}
export default compose(withApi(false, false), withStyles(UploadersList.styles))(UploadersList)
