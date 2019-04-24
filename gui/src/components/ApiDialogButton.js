import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, IconButton, Dialog, DialogTitle, DialogContent, DialogActions, Button } from '@material-ui/core'
import CodeIcon from '@material-ui/icons/Code'
import ReactJson from 'react-json-view'

class ApiDialogUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    data: PropTypes.any.isRequired,
    title: PropTypes.string,
    onClose: PropTypes.func
  }

  static styles = (theme) => ({
    content: {
      paddingBottom: 0
    },
    raw: {
      margin: 0, padding: 0
    }
  })

  state = {
    showRaw: false
  }

  constructor(props) {
    super(props)
    this.handleToggleRaw = this.handleToggleRaw.bind(this)
  }

  handleToggleRaw() {
    this.setState({showRaw: !this.state.showRaw})
  }

  render() {
    const { classes, title, data, onClose, ...dialogProps } = this.props
    const { showRaw } = this.state

    return (
      <Dialog {...dialogProps}>
        <DialogTitle>{title || 'API'}</DialogTitle>
        <DialogContent classes={{root: classes.content}}>
          {showRaw
            ? <code>
              <pre className={classes.raw}>
                {JSON.stringify(data, null, 4)}
              </pre>
            </code> : <ReactJson
              src={data}
              enableClipboard={false}
              collapsed={2}
              displayObjectSize={false}
            />
          }
        </DialogContent>
        <DialogActions>
          <Button onClick={this.handleToggleRaw} color="primary">
            {showRaw ? 'show tree' : 'show raw JSON'}
          </Button>
          <Button onClick={onClose} color="primary">
            Close
          </Button>
        </DialogActions>
      </Dialog>
    )
  }
}

export const ApiDialog = withStyles(ApiDialogUnstyled.styles)(ApiDialogUnstyled)

class ApiDialogButtonUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    data: PropTypes.any.isRequired,
    title: PropTypes.string
  }

  static styles = theme => ({
    root: {}
  })

  state = {
    showDialog: false
  }

  render() {
    const { classes, ...dialogProps } = this.props
    const { showDialog } = this.state

    return (
      <div className={classes.root}>
        <IconButton onClick={() => this.setState({showDialog: true})}>
          <CodeIcon />
        </IconButton>
        <ApiDialog
          {...dialogProps} open={showDialog}
          onClose={() => this.setState({showDialog: false})}
        />
      </div>
    )
  }
}

export default withStyles(ApiDialogButtonUnstyled.styles)(ApiDialogButtonUnstyled)
