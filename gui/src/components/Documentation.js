import React, { Component } from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core';
import { apiBase } from '../config';


class Documentation extends Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
  }

  static styles = theme => ({
    root: {
      position: 'relative',
    },
    content: {
      position: 'absolute',
      left: -theme.spacing.unit * 3,
      top: -theme.spacing.unit * 3
    }
  })

  render() {
    const {classes} = this.props
    return (
      <div className={classes.root}>
        <div className={classes.content}>
          <iframe title="documentation"
            frameBorder={0} width="700" height={window.innerHeight - 64}
            src={`${apiBase}/docs/index.html`}
          />
        </div>
      </div>
    )
  }
}

export default withStyles(Documentation.styles)(Documentation)
