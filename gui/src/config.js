import repo from '@material-ui/core/colors/deepPurple';
import archive from '@material-ui/core/colors/teal';
import enc from '@material-ui/core/colors/amber';
import secondary from '@material-ui/core/colors/blueGrey';
import { createMuiTheme } from '@material-ui/core';

export const apiBase = 'http://localhost:5000'

export const repoTheme = createMuiTheme({
  palette: {
    primary: repo,
    secondary: secondary,
  }
});

export const archiveTheme = createMuiTheme({
  palette: {
    primary: archive,
    secondary: secondary,
  }
});

export const encTheme = createMuiTheme({
  palette: {
    primary: enc,
    secondary: secondary,
  }
});