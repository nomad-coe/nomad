import React from 'react';
import { MuiThemeProvider } from '@material-ui/core/styles';
import { repoTheme } from '../config';
import Navigation from './Navigation';
import { BrowserRouter, Switch, Route } from 'react-router-dom';
import Uploads from './Uploads'
import ArchiveCalc from './ArchiveCalc';

function App() {
  return (
    <MuiThemeProvider theme={repoTheme}>
      <BrowserRouter>
        <Navigation>
          <Switch>
            <Route exact path="/" render={() => <div>Home</div>} />
            <Route path="/repo" render={() => <div>Browse</div>} />
            <Route path="/upload" component={Uploads} />
            <Route exact path="/archive" render={() => <div>Archive</div>} />
            <Route path="/archive/:uploadHash/:calcHash" component={ArchiveCalc} />
            <Route path="/profile" render={() => <div>Profile</div>} />
            <Route path="/documentation" render={() => <div>Docs</div>} />
            <Route render={() => <div>Not found</div>}  />
          </Switch>
        </Navigation>
      </BrowserRouter>
    </MuiThemeProvider>
  );
}

export default App;
