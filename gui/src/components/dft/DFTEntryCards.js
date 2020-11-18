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
import { withStyles, Card, CardHeader, CardContent } from '@material-ui/core'
import RawFiles from '../entry/RawFiles'

class DFTEntryCards extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired
  }

  static styles = theme => ({
    root: {}
  })

  render() {
    const { classes, ...props } = this.props

    return (
      <Card className={classes.root}>
        <CardHeader title="Raw files" />
        <CardContent classes={{root: classes.cardContent}}>
          <RawFiles {...props} />
        </CardContent>
      </Card>
    )
  }
}

export default withStyles(DFTEntryCards.styles)(DFTEntryCards)
