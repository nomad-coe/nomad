import React from 'react'
import { withStyles, Button, Collapse } from '@material-ui/core'
import Markdown from './Markdown'
import PropTypes, { instanceOf } from 'prop-types'
import { Cookies, withCookies } from 'react-cookie'
import classNames from 'classnames'

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

export class Help extends React.Component {
  static propTypes = {
    children: PropTypes.any,
    cookie: PropTypes.string.isRequired
  }

  render() {
    const { children, cookie } = this.props

    return (
      <HelpContext.Consumer>{
        help => (
          <Collapse in={help.isOpen(cookie)}>
            <GotIt onGotIt={() => help.gotIt(cookie)}>
              {children}
            </GotIt>
          </Collapse>
        )
      }</HelpContext.Consumer>
    )
  }
}

class GotItUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    onGotIt: PropTypes.func.isRequired,
    children: PropTypes.node.isRequired,
    color: PropTypes.oneOf(['primary', 'error']).isRequired
  }

  static defaultProps = {
    color: 'primary'
  }

  static styles = theme => ({
    root: {
      marginTop: theme.spacing.unit * 2,
      marginBottom: theme.spacing.unit * 2,
      borderRadius: theme.spacing.unit * 0.5,
      display: 'flex',
      flexDirection: 'row',
      alignItems: 'center'
    },
    rootPrimary: {
      border: `1px solid ${theme.palette.primary.main}`
    },
    rootError: {
      border: `1px solid ${theme.palette.error.main}`
    },
    content: {
      padding: theme.spacing.unit * 2,
      flex: '1 1 auto'
    },
    actions: {
      padding: theme.spacing.unit,
      paddingLef: 0,
      flex: '0 0 auto'
    }
  })

  render() {
    const { classes, children, onGotIt, color } = this.props
    const rootClassName = classNames(classes.root, {
      [classes.rootPrimary]: color === 'primary',
      [classes.rootError]: color === 'error'
    })
    return (
      <div className={rootClassName}>
        <div className={classes.content}>
          <Markdown>
            {children}
          </Markdown>
        </div>
        <div className={classes.actions}>
          <Button color="primary" onClick={onGotIt}>Got it</Button>
        </div>
      </div>
    )
  }
}

const GotIt = withStyles(GotItUnstyled.styles)(GotItUnstyled)

export const HelpProvider = withCookies(HelpProviderComponent)
