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
import { withStyles, Typography, Link } from '@material-ui/core'
import { compose } from 'recompose'
import { withApi, DoesNotExist } from '../api'
import { withRouter } from 'react-router'
import qs from 'qs'

class EntryQuery extends React.Component {
  static styles = theme => ({
    root: {
      padding: theme.spacing(3)
    }
  })

  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    history: PropTypes.object.isRequired,
    location: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired
  }

  static defaultState = {
    doesNotExist: false,
    queryParams: null
  }

  state = {...EntryQuery.defaultState}

  update() {
    const { api, history, location } = this.props

    let queryParams = null
    if (location && location.search) {
      queryParams = qs.parse(location.search.substring(1))
    }
    api.search({...queryParams}).then(data => {
      if (data.results && data.results.length > 0) {
        const { calc_id, upload_id } = data.results[0]
        history.push(`/entry/id/${upload_id}/${calc_id}`)
      } else {
        throw new DoesNotExist()
      }
    }).catch(error => {
      if (error.name === 'DoesNotExist') {
        this.setState({doesNotExist: true, queryParams: queryParams})
      } else {
        this.props.raiseError(error)
      }
    })
  }

  componentDidMount() {
    this.update()
  }

  componentDidUpdate(prevProps) {
    if (prevProps.location.search !== this.props.location.search || prevProps.api !== this.props.api) {
      this.setState({...EntryQuery.defaultState})
      this.update()
    }
  }

  render() {
    const { classes, api } = this.props
    const { doesNotExist, queryParams } = this.state

    let message = 'loading ...'

    if (doesNotExist) {
      console.log(queryParams)
      if (queryParams && queryParams['external_id'] && queryParams['external_id'].startsWith('mp-')) {
        message = <React.Fragment>
          This particular calculation <Link href={`https://materialsproject.org/tasks/${queryParams['external_id']}#`}>
            {queryParams['external_id']}
          </Link> has not yet been provided to NOMAD by the Materials Project.
        </React.Fragment>
      } else if (api.isLoggedIn) {
        message = `
            This URL points to an entry that either does not exist, or that you are not
            authorized to see.`
      } else {
        message = `
            This URL points to an entry that either does not exist, or that is not
            publically visibile. Please login; you might be authorized to view it.`
      }
    }

    return (
      <Typography className={classes.root}>{message}</Typography>
    )
  }
}

export default compose(withRouter, withApi(false), withStyles(EntryQuery.styles))(EntryQuery)
