import React from 'react'
import PropTypes from 'prop-types'
import { Typography, withStyles } from '@material-ui/core'

class UploadCommand extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    uploadCommand: PropTypes.string.isRequired
  }

  static styles = theme => ({
    root: {
      margin: theme.spacing.unit * 2
    },
    uploadCommand: {
      fontFamily: '\'Roboto mono\', monospace',
      marginTop: theme.spacing.unit * 2
    }
  })

  render() {
    const { classes, uploadCommand } = this.props
    return (
      <div className={classes.root}>
        <Typography>Copy and use the following command. Don't forget to replace the file name.:</Typography>
        <Typography className={classes.uploadCommand}>{uploadCommand}</Typography>
      </div>
    )
  }
}

export default withStyles(UploadCommand.styles)(UploadCommand)
