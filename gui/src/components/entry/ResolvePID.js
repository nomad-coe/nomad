import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Typography } from '@material-ui/core'
import { compose } from 'recompose'
import { withApi } from '../api'
import { withRouter } from 'react-router'

class ResolvePID extends React.Component {
  static styles = theme => ({
    root: {
      padding: theme.spacing.unit * 3
    }
  })

  static propTypes = {
    classes: PropTypes.object.isRequired,
    pid: PropTypes.string.isRequired,
    api: PropTypes.object.isRequired,
    history: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired
  }

  static defaultState = {
    doesNotExist: false
  }

  state = {...ResolvePID.defaultState}

  update() {
    const { pid, api, history } = this.props
    api.resolvePid(pid).then(entry => {
      history.push(`/entry/id/${entry.upload_id}/${entry.calc_id}`)
    }).catch(error => {
      this.setState({calcData: null})
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
    if (prevProps.pid !== this.props.pid || prevProps.api !== this.props.api) {
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
