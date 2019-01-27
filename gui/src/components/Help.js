import React from 'react'
import { withStyles, Button } from '@material-ui/core'
import Markdown from './Markdown'
import PropTypes, { instanceOf } from 'prop-types'
import { Cookies, withCookies } from 'react-cookie'

export const HelpContext = React.createContext()

class HelpProviderComponent extends React.Component {
  static propTypes = {
    children: PropTypes.oneOfType([
      PropTypes.arrayOf(PropTypes.node),
      PropTypes.node
    ]).isRequired,
    cookies: instanceOf(Cookies).isRequired
  }

  state = {
    helpCookies: [],
    allHelpCookies: [],
    allClosed: () => this.state.helpCookies.length === this.state.allHelpCookies.length,
    someClosed: () => this.state.helpCookies.length !== 0,
    isOpen: (cookie) => {
      if (this.state.allHelpCookies.indexOf(cookie) === -1) {
        this.state.allHelpCookies.push(cookie)
      }
      return this.state.helpCookies.indexOf(cookie) === -1
    },
    gotIt: (cookie) => {
      const updatedHelpCookies = [...this.state.helpCookies, cookie]
      this.props.cookies.set('help', updatedHelpCookies)
      this.setState({helpCookies: updatedHelpCookies})
    },
    switchHelp: () => {
      const updatedCookies = this.state.someClosed() ? [] : this.state.allHelpCookies
      this.setState({helpCookies: updatedCookies})
      this.props.cookies.set('help', updatedCookies)
    }
  }

  componentDidMount() {
    this.setState({helpCookies: this.props.cookies.get('help') || []})
  }

  render() {
    return (
      <HelpContext.Provider value={this.state}>
        {this.props.children}
      </HelpContext.Provider>
    )
  }
}

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

export const HelpProvider = withCookies(HelpProviderComponent)
export const Help = withStyles(HelpComponent.styles)(HelpComponent)
