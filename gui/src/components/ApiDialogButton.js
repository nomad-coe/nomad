import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, IconButton, Dialog, DialogTitle, DialogContent, DialogActions, Button } from '@material-ui/core'
import CodeIcon from '@material-ui/icons/Code'
import ReactJson from 'react-json-view'

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
    const { classes, title, data } = this.props
    const { showDialog } = this.state

    return (
      <div className={classes.root}>
        <IconButton onClick={() => this.setState({showDialog: true})}>
          <CodeIcon />
        </IconButton>
        <Dialog open={showDialog}>
          <DialogTitle>{title || 'API'}</DialogTitle>
          <DialogContent>
            <ReactJson
              src={data}
              enableClipboard={false}
              collapsed={2}
              displayObjectSize={false}
            />
          </DialogContent>
          <DialogActions>
            <Button onClick={() => this.setState({showDialog: false})} color="primary">
                Close
            </Button>
          </DialogActions>
        </Dialog>
      </div>
    )
  }
}

export default withStyles(ApiDialogButtonUnstyled.styles)(ApiDialogButtonUnstyled)
