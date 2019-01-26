import React from 'react'
import { withStyles, Button } from '@material-ui/core'
import Markdown from './Markdown'
import PropTypes from 'prop-types'

export class HelpManager {
  isOpen(helpKey) {
    return true
  }
  gotIt(helpKey) {

  }
}

export const HelpContext = React.createContext(new HelpManager())

class HelpComponent extends React.Component {
  static styles = theme => ({
    root: {
      marginTop: theme.spacing.unit * 2,
      marginBottom: theme.spacing.unit * 2,
      borderRadius: theme.spacing.unit * 0.5,
      border: `1px solid ${theme.palette.primary.main}`,
      display: 'flex',
      flexDirection: 'row',
      alignItems: 'center'
    },
    content: {
      paddingLeft: theme.spacing.unit * 2,
      flex: '1 1 auto'
    },
    actions: {
      padding: theme.spacing.unit,
      flex: '0 0 auto'
    }
  })

  static propTypes = {
    classes: PropTypes.object.isRequired,
    children: PropTypes.any,
    cookie: PropTypes.string.isRequired
  }

  render() {
    const { classes, children, cookie } = this.props

    return (
      <HelpContext.Consumer>{
        help => (
          help.isOpen(cookie)
            ? <div className={classes.root}>
              <div className={classes.content}>
                <Markdown>
                  {children}
                </Markdown>
              </div>
              <div className={classes.actions}>
                <Button color="primary" onClick={() => help.gotIt(cookie)}>Got it</Button>
              </div>
            </div> : ''
        )
      }</HelpContext.Consumer>
    )
  }
}

export const Help = withStyles(HelpComponent.styles)(HelpComponent)
