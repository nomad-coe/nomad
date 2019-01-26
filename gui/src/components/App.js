import React from 'react'
import { MuiThemeProvider } from '@material-ui/core/styles'
import { genTheme, appBase } from '../config'
import Navigation from './Navigation'
import { BrowserRouter, Switch, Route } from 'react-router-dom'
import Uploads from './Uploads'
import Repo from './Repo'
import Documentation from './Documentation'
import Development from './Development'
import Home from './Home'
import { HelpContext } from './Help'
import { withCookies, Cookies } from 'react-cookie'
import { instanceOf } from 'prop-types'

class App extends React.Component {
  static propTypes = {
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
      <MuiThemeProvider theme={genTheme}>
        <BrowserRouter basename={appBase}>
          <HelpContext.Provider value={this.state}>
            <Navigation>
              <Switch>
                <Route exact path="/" component={Home} />
                <Route exact path="/repo" component={Repo} />
                {/* <Route path="/repo/:uploadId/:calcId" component={RepoCalc} /> */}
                <Route path="/upload" component={Uploads} />
                <Route exact path="/archive" render={() => <div>Archive</div>} />
                {/* <Route path="/archive/:uploadId/:calcId" component={ArchiveCalc} /> */}
                <Route path="/enc" render={() => <div>{'In the future, you\'ll see charts\'n\'stuff for your calculations and materials.'}</div>} />
                <Route path="/analytics" render={() => <div>{'In the future, you\'ll see analytics notebooks here.'}</div>} />
                <Route path="/profile" render={() => <div>Profile</div>} />
                <Route path="/docs" component={Documentation} />
                <Route path="/dev" component={Development} />
                <Route render={() => <div>Not found</div>} />
              </Switch>
            </Navigation>
          </HelpContext.Provider>
        </BrowserRouter>
      </MuiThemeProvider>
    )
  }
}

export default withCookies(App)
