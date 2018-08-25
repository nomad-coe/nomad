import React from 'react';
import PropTypes from 'prop-types';
import { MuiThemeProvider, IconButton, withStyles } from "@material-ui/core";
import RepoIcon from '@material-ui/icons/Cloud';
import ArchiveIcon from '@material-ui/icons/Storage';
import EncIcon from '@material-ui/icons/Assessment';
import { repoTheme, archiveTheme, encTheme } from '../config';
import Link from 'react-router-dom/Link';


class CalcLink extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    uploadHash: PropTypes.string.isRequired,
    calcHash: PropTypes.string.isRequired,
  }
  static styles = theme => ({
      root: {
        overflow: 'hidden',
        whiteSpace: 'nowrap',
        textAlign: 'right'
      }
  });

  render() {
    const { uploadHash, calcHash, classes } = this.props
    const archiveId = `${uploadHash}/${calcHash}`

    return (
      <div className={classes.root}>
        <MuiThemeProvider theme={repoTheme}>
          <IconButton color="primary" component={Link} to={`/repo/${archiveId}`}><RepoIcon /></IconButton>
        </MuiThemeProvider>
        <MuiThemeProvider theme={archiveTheme}>
          <IconButton color="primary" component={Link} to={`/archive/${archiveId}`}><ArchiveIcon /></IconButton>
        </MuiThemeProvider>
        <MuiThemeProvider theme={encTheme}>
          <IconButton color="primary" component={Link} to={`/enc/${archiveId}`}><EncIcon /></IconButton>
        </MuiThemeProvider>
      </div>
    )
  }
}

export default withStyles(CalcLink.styles)(CalcLink);