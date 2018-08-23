import React from 'react';
import { MuiThemeProvider, createMuiTheme } from '@material-ui/core/styles';
import deepPurple from '@material-ui/core/colors/deepPurple';
import orange from '@material-ui/core/colors/orange';
import Navigation from './Navigation';
import { BrowserRouter, Switch, Route } from 'react-router-dom';
import Uploads from './Uploads'

const theme = createMuiTheme({
  palette: {
    primary: deepPurple,
    secondary: orange,
  },
});

function App() {
  return (
    <MuiThemeProvider theme={theme}>
      <BrowserRouter>
        <Navigation>
          <Switch>
            <Route exact path="/" render={() => <div>Home</div>} />
            <Route path="/browse" render={() => <div>Browse</div>} />
            <Route path="/upload" component={Uploads} />
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
