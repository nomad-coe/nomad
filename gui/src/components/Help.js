import React from 'react'
import { withStyles, Button, IconButton, Dialog, DialogTitle, DialogContent, DialogActions, Tooltip } from '@material-ui/core'
import Markdown from './Markdown'
import PropTypes from 'prop-types'
import HelpIcon from '@material-ui/icons/Help'
import { compose } from 'recompose'
import { withDomain } from './domains'

export const HelpContext = React.createContext()

class HelpDialogUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    title: PropTypes.string,
    content: PropTypes.func.isRequired,
    icon: PropTypes.node,
    maxWidth: PropTypes.string,
    domains: PropTypes.object.isRequired
  }

  static styles = theme => ({
    root: {}
  })

  state = {
    isOpen: false
  }

  constructor(props) {
    super(props)
    this.handleOpen = this.handleOpen.bind(this)
    this.handleClose = this.handleClose.bind(this)
  }

  handleClose() {
    this.setState({isOpen: false})
  }

  handleOpen() {
    this.setState({isOpen: true})
  }

  render() {
    const {classes, title, content, icon, maxWidth, domains, ...rest} = this.props
    return (
      <div className={classes.root}>
        <Tooltip title={title}>
          <IconButton {...rest} onClick={this.handleOpen}>
            {icon || <HelpIcon/>}
          </IconButton>
        </Tooltip>
        <Dialog
          maxWidth={maxWidth}
          onClose={this.handleClose}
          open={this.state.isOpen}
        >
          <DialogTitle>{title || 'Help'}</DialogTitle>
          <DialogContent>
            <Markdown>{content(domains)}</Markdown>
          </DialogContent>
          <DialogActions>
            <Button onClick={() => this.handleClose()} color="primary">
              Close
            </Button>
          </DialogActions>
        </Dialog>
      </div>
    )
  }
}

export default compose(withDomain, withStyles(HelpDialogUnstyled.styles))(HelpDialogUnstyled)
