import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Typography } from '@material-ui/core'
import { compose } from 'recompose'
import { withApi } from '../api'
import { matchPath } from 'react-router'

class ResolveDOI extends React.Component {
  static styles = theme => ({
    root: {
      padding: theme.spacing.unit * 3
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
