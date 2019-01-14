import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import Paper from '@material-ui/core/Paper'
import api from '../api'
import { Popover } from '@material-ui/core'


class CalcProcLogPopper extends React.Component {

  static propTypes = {
    classes: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    uploadId: PropTypes.string.isRequired,
    calcId: PropTypes.string.isRequired,
    open: PropTypes.bool,
    onClose: PropTypes.func,
    anchorEl: PropTypes.any
  }

  static styles = theme => ({
    paper: {
      padding: theme.spacing.unit * 2,
    },
  })

  state = {
    logs: null
  }

  componentDidMount() {
    const {uploadId, calcId} = this.props
    api.calcProcLog(uploadId, calcId).then(logs => {
      if (logs && logs !== '') {
        this.setState({logs: logs})
      }
    }).catch(error => {
      this.setState({data: null})
      this.props.raiseError(error)
    })
  }

  render() {
    const { classes, open, anchorEl, onClose } = this.props
    const { logs } = this.state
    return (
      <div>
        <Popover
          open={open}
          anchorEl={anchorEl}
          onClose={onClose}
          anchorOrigin={{
            vertical: 'center',
            horizontal: 'center',
          }}
          transformOrigin={{
            vertical: 'center',
            horizontal: 'center',
          }}
        >
          <Paper className={classes.paper}>
            <pre>
              {logs ? logs : 'loading...'}
            </pre>
          </Paper>
        </Popover>
      </div>
    )
  }
}

export default withStyles(CalcProcLogPopper.styles)(CalcProcLogPopper)