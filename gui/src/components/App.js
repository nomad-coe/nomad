import React from 'react'
import { MuiThemeProvider } from '@material-ui/core/styles'
import { genTheme, appBase } from '../config'
import Navigation from './Navigation'
import { BrowserRouter, Switch, Route } from 'react-router-dom'
import Uploads from './Uploads'
import ArchiveCalc from './ArchiveCalc'
import RepoCalc from './RepoCalc'
import Repo from './Repo'
import Documentation from './Documentation'
import Development from './Development'
import Home from './Home'

function App() {
  return (
    <MuiThemeProvider theme={genTheme}>
      <BrowserRouter basename={appBase}>
        <Navigation>
          <Switch>
            <Route exact path="/" component={Home} />
            <Route exact path="/repo" component={Repo} />
            <Route path="/repo/:uploadId/:calcId" component={RepoCalc} />
            <Route path="/upload" component={Uploads} />
            <Route exact path="/archive" render={() => <div>Archive</div>} />
            <Route path="/archive/:uploadId/:calcId" component={ArchiveCalc} />
            <Route path="/enc" render={() => <div>{'In the future, you\'ll see charts\'n\'stuff for your calculations and materials.'}</div>} />
            <Route path="/analytics" render={() => <div>{'In the future, you\'ll see analytics notebooks here.'}</div>} />
            <Route path="/profile" render={() => <div>Profile</div>} />
            <Route path="/docs" component={Documentation} />
            <Route path="/dev" component={Development} />
            <Route render={() => <div>Not found</div>} />
          </Switch>
        </Navigation>
      </BrowserRouter>
    </MuiThemeProvider>
  )
}

export default App
