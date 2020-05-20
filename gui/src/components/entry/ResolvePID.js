import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Typography } from '@material-ui/core'
import { compose } from 'recompose'
import { withApi } from '../api'
import { withRouter, matchPath } from 'react-router'

class ResolvePID extends React.Component {
  static styles = theme => ({
    root: {
      padding: theme.spacing(3)
    }
  })

  static propTypes = {
    classes: PropTypes.object.isRequired,
    location: PropTypes.object.isRequired,
    match: PropTypes.object.isRequired,
    history: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired
  }

  static defaultState = {
    doesNotExist: false
  }

  state = {...ResolvePID.defaultState}

  update() {
    const { location, match, api, history } = this.props
    const pidMatch = matchPath(location.pathname, {
      path: `${match.path}/:pid/:handle?`
    })
    let { pid, handle } = pidMatch.params
    pid = handle ? pid + '/' + handle : pid

    api.resolvePid(pid).then(entry => {
      history.push(`/entry/id/${entry.upload_id}/${entry.calc_id}`)
    }).catch(error => {
      if (error.name === 'DoesNotExist') {
        this.setState({doesNotExist: true})
      } else {
        this.props.raiseError(error)
      }
    })
  }

  componentDidMount() {
    this.update()
  }

  componentDidUpdate(prevProps) {
    if (prevProps.location.pathname !== this.props.location.pathname || prevProps.api !== this.props.api) {
      this.setState({...ResolvePID.defaultState})
      this.update()
    }
  }

  render() {
    const { classes, api } = this.props
    const { doesNotExist } = this.state

    let message = 'loading ...'

    if (doesNotExist) {
      if (api.isLoggedIn) {
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

export default compose(withRouter, withApi(false), withStyles(ResolvePID.styles))(ResolvePID)
