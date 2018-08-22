import React from 'react';
import { MuiThemeProvider, createMuiTheme } from '@material-ui/core/styles';
import blue from '@material-ui/core/colors/blue';
import orange from '@material-ui/core/colors/orange';
import Navigation from './Navigation';
import { BrowserRouter, Switch, Route } from 'react-router-dom';

const theme = createMuiTheme({
  palette: {
    primary: blue,
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
            <Route path="/upload" render={() => <div>Upload</div>} />
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
