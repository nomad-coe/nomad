import React from 'react'
import { MuiThemeProvider } from '@material-ui/core/styles'
import { genTheme, appBase } from '../config'
import Navigation from './Navigation'
import { BrowserRouter, Switch, Route } from 'react-router-dom'
import Uploads from './uploads/Uploads'
import SearchPage from './search/SearchPage'
import Development from './Development'
import Home from './Home'
import { HelpProvider } from './help'
import { ApiProvider } from './api'
import { ErrorSnacks } from './errors'
import Calc from './entry/Calc'

export default class App extends React.Component {
  constructor(props) {
    super(props)
    this.renderChildren.bind(this)
  }

  routes = {
    'home': {
      exact: true,
      path: '/',
      render: props => <Home {...props} />
    },
    'search': {
      exact: true,
      path: '/search',
      render: props => <SearchPage {...props} />
    },
    'searchEntry': {
      path: '/search/:uploadId/:calcId',
      render: props => {
        const { match, ...rest } = props
        if (match && match.params.uploadId && match.params.calcId) {
          return (<Calc {...rest} uploadId={match.params.uploadId} calcId={match.params.calcId} />)
        } else {
          return ''
        }
      }
    },
    'uploads': {
      exact: true,
      path: '/uploads',
      render: props => <Uploads {...props} />
    },
    'uploadedEntry': {
      path: '/uploads/:uploadId/:calcId',
      render: props => {
        const { match, ...rest } = props
        if (match && match.params.uploadId && match.params.calcId) {
          return (<Calc {...rest} uploadId={match.params.uploadId} calcId={match.params.calcId} />)
        } else {
          return ''
        }
      }
    },
    'dev': {
      exact: true,
      path: '/dev',
      render: props => <Development {...props} />
    }
  }

  renderChildren(routeKey, props) {
    // const { match, ...rest } = props

    return (
      <div>
        {Object.keys(this.routes).map(route => (
          <div key={route} style={{display: routeKey === route ? 'block' : 'none'}}>
            {this.routes[route].render(props)}
          </div>
        ))}
      </div>
    )
  }

  render() {
    return (
      <MuiThemeProvider theme={genTheme}>
        <ErrorSnacks>
          <BrowserRouter basename={appBase}>
            <HelpProvider>
              <ApiProvider>
                <Navigation>
                  <Switch>
                    {Object.keys(this.routes).map(route => (
                      // eslint-disable-next-line react/jsx-key
                      <Route key={'nop'}
                        // eslint-disable-next-line react/no-children-prop
                        children={props => this.renderChildren(route, props)}
                        exact={this.routes[route].exact}
                        path={this.routes[route].path} />
                    ))}
                  </Switch>
                </Navigation>
              </ApiProvider>
            </HelpProvider>
          </BrowserRouter>
        </ErrorSnacks>
      </MuiThemeProvider>
    )
  }
}
