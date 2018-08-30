import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import Markdown from './Markdown'
import gitInfo from '../gitinfo'

class Development extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired
  }

  static styles = theme => ({
    root: {},
  })

  render() {
    const { classes } = this.props

    return (
      <div className={classes.root}>
        <Markdown>{`
          # Build info
          - version: \`${gitInfo.version}\`
          - ref: \`${gitInfo.ref}\`
          - last commit message: *${gitInfo.log}*
        `}</Markdown>
      </div>
    )
  }
}

export default withStyles(Development.styles)(Development)
