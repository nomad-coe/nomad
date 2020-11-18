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
import { withStyles, Typography } from '@material-ui/core'
import { compose } from 'recompose'
import { withApi } from '../api'
import { matchPath } from 'react-router'

class ResolveDOI extends React.Component {
  static styles = theme => ({
    root: {
      padding: theme.spacing(3)
    }
  })

  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    location: PropTypes.object.isRequired,
    match: PropTypes.object.isRequired,
    history: PropTypes.object.isRequired
  }

  update() {
    const { location, match, api, history, raiseError } = this.props
    const doiMatch = matchPath(location.pathname, {
      path: `${match.path}/:doi*`
    })
    let { doi } = doiMatch.params

    api.resolveDoi(doi).then(dataset => {
      history.push(`/dataset/id/${dataset.dataset_id}`)
    }).catch(raiseError)
  }

  componentDidMount() {
    this.update()
  }

  componentDidUpdate(prevProps) {
    if (prevProps.location.pathname !== this.props.location.pathname || prevProps.api !== this.props.api) {
      this.update()
    }
  }

  render() {
    const { classes } = this.props

    return (
      <Typography className={classes.root}>loading ...</Typography>
    )
  }
}

export default compose(withApi(false), withStyles(ResolveDOI.styles))(ResolveDOI)
