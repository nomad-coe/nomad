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
import Markdown from '../Markdown'

class EMSEntryRawView extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    data: PropTypes.object.isRequired,
    loading: PropTypes.bool
  }

  static styles = theme => ({
    description: {
      marginBottom: theme.spacing(3)
    }
  })

  render() {
    const { classes, data, ...props } = this.props

    return (
      <Card>
        <CardHeader title="Raw Data and Meta Data Files" />
        <CardContent classes={{root: classes.cardContent}}>
          <Markdown classes={{root: classes.description}}>{`
            The data for this experiment is externally stored and managed. Download the raw experiment data:
            [${data.ems && data.ems.repository_url}](${data.ems && data.ems.entry_repository_url}).

            The meta data describing this experiment in its original format, can be
            downloaded here directly:
          `}</Markdown>
          <RawFiles data={data} {...props} />
        </CardContent>
      </Card>
    )
  }
}

export default withStyles(EMSEntryRawView.styles)(EMSEntryRawView)
