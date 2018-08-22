import React from 'react';
import PropTypes from 'prop-types';
import { withStyles } from '@material-ui/core/styles';
import Drawer from '@material-ui/core/Drawer';
import AppBar from '@material-ui/core/AppBar';
import Toolbar from '@material-ui/core/Toolbar';
import Typography from '@material-ui/core/Typography';
import Divider from '@material-ui/core/Divider';
import MenuList from '@material-ui/core/MenuList';
import MenuItem from '@material-ui/core/MenuItem';
import ListItemIcon from '@material-ui/core/ListItemIcon';
import ListItemText from '@material-ui/core/ListItemText';
import BackupIcon from '@material-ui/icons/Backup';
import SearchIcon from '@material-ui/icons/Search';
import AccountIcon from '@material-ui/icons/AccountCircle';
import DocumentationIcon from '@material-ui/icons/Help';
import HomeIcon from '@material-ui/icons/Home';
import { Link, withRouter } from 'react-router-dom';
import { compose } from 'recompose'
import { Avatar } from '@material-ui/core';

const drawerWidth = 200;

const toolbarTitles = {
    '/': 'Welcome',
    '/browse': 'Search, View, and Download Data',
    '/upload': 'Upload Your Own Data',
    '/profile': 'Your Profile',
    '/documentation': 'Documentation'
}

const styles = theme => ({
  root: {
    flexGrow: 1
  },
  flex: {
    flexGrow: 1,
  },
  appFrame: {
    height: 440,
    zIndex: 1,
    overflow: 'hidden',
    position: 'relative',
    display: 'flex',
    width: '100%',
  },
  appBar: {
    width: `calc(100% - ${drawerWidth}px)`,
    marginLeft: drawerWidth,
  },
  drawerPaper: {
    position: 'relative',
    width: drawerWidth,
  },
  content: {
    flexGrow: 1,
    backgroundColor: theme.palette.background.default,
    padding: theme.spacing.unit * 3,
  },
  toolbar: theme.mixins.toolbar,
  link: {
      textDecoration: 'none',
      color: theme.palette.text.primary
  }
});

function ClippedDrawer(props) {
  const { classes, children, location: { pathname } } = props;

  const title = (
    <Toolbar>
      <Typography variant="title" color="inherit" noWrap>
        nomad xt
      </Typography>
    </Toolbar>
  )

  const menu = (
    <div>
      <MenuList>
        <MenuItem component={Link} to="/" selected={ '/' === pathname }>
          <ListItemIcon>
            <HomeIcon />
          </ListItemIcon>
          <ListItemText inset primary="Home"/>
        </MenuItem>
        <MenuItem component={Link} to="/browse" selected={ '/browse' === pathname }>
          <ListItemIcon>
            <SearchIcon />
          </ListItemIcon>
          <ListItemText inset primary="Browse"/>
        </MenuItem>
        <MenuItem component={Link} to="/upload" selected={ '/upload' === pathname }>
          <ListItemIcon>
            <BackupIcon />
          </ListItemIcon>
          <ListItemText inset primary="Upload"/>
        </MenuItem>
      </MenuList>
      <Divider/>
      <MenuList>
        <MenuItem component={Link} to="/profile" selected={ '/profile' === pathname }>
          <ListItemIcon>
            <AccountIcon />
          </ListItemIcon>
          <ListItemText inset primary="Profil"/>
        </MenuItem>
        <MenuItem component={Link} to="/documentation" selected={ '/documentation' === pathname }>
          <ListItemIcon>
            <DocumentationIcon />
          </ListItemIcon>
          <ListItemText inset primary="Documentation"/>
        </MenuItem>
      </MenuList>
    </div>
  )

  const drawer = (
    <Drawer variant="permanent" classes={{ paper: classes.drawerPaper, }} anchor="left">
      <div className={classes.toolbar}>
        {title}
      </div>
      <Divider />
      {menu}
    </Drawer>
  );

  const app = (
    <div className={classes.root}>
      <div className={classes.appFrame}>
        <AppBar
          position="absolute"
          className={classes.appBar}
        >
          <Toolbar>
            <Typography variant="title" color="inherit" noWrap className={classes.flex}>
              {toolbarTitles[pathname]}
            </Typography>
            <Avatar src='/me.jpg'/>
          </Toolbar>
        </AppBar>
        {drawer}
        <main className={classes.content}>
          <div className={classes.toolbar} />
          {children}
        </main>
      </div>
    </div>
  );

  return app;
}

ClippedDrawer.propTypes = {
  classes: PropTypes.object.isRequired,
};

export default compose(withRouter, withStyles(styles))(ClippedDrawer)

