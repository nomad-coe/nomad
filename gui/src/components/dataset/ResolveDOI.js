import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Typography } from '@material-ui/core'
import { compose } from 'recompose'
import { withApi } from '../api'
import { withRouter } from 'react-router'

class ResolveDOI extends React.Component {
  static styles = theme => ({
    root: {
      padding: theme.spacing.unit * 3
    }
  })

  static propTypes = {
    classes: PropTypes.object.isRequired,
    doi: PropTypes.string.isRequired,
    api: PropTypes.object.isRequired,
    history: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired
  }

  update() {
    const { doi, api, history, raiseError } = this.props
    api.resolveDoi(doi).then(dataset => {
      history.push(`/dataset/id/${dataset.dataset_id}`)
    }).catch(raiseError)
  }

  componentDidMount() {
    this.update()
  }

  componentDidUpdate(prevProps) {
    if (prevProps.doi !== this.props.doi || prevProps.api !== this.props.api) {
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

export default compose(withRouter, withApi(false), withStyles(ResolveDOI.styles))(ResolveDOI)
