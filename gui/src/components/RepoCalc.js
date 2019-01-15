import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Paper } from '@material-ui/core'
import Markdown from './Markdown'
import { withErrors } from './errors'
import { compose } from 'recompose'
import RepoCalcView from './RepoCalcView'

class RepoCalc extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    match: PropTypes.object.isRequired
  }

  static styles = theme => ({
    root: {},
    calcData: {
      padding: theme.spacing.unit
    }
  })

  render() {
    const { classes, match, raiseError } = this.props
    const {uploadId, calcId} = match.params

    return (
      <div className={classes.root}>
        <Markdown>{`
          ## The Repository – Raw Code Data
        `}</Markdown>
        <Paper className={classes.calcData}>
          <RepoCalcView uploadId={uploadId} calcId={calcId} raiseError={raiseError} />
        </Paper>
      </div>

    )
  }
}

export default compose(withErrors, withStyles(RepoCalc.styles))(RepoCalc)
