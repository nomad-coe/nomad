import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import Drawer from '@material-ui/core/Drawer'
import AppBar from '@material-ui/core/AppBar'
import Toolbar from '@material-ui/core/Toolbar'
import Typography from '@material-ui/core/Typography'
import Divider from '@material-ui/core/Divider'
import MenuList from '@material-ui/core/MenuList'
import MenuItem from '@material-ui/core/MenuItem'
import ListItemIcon from '@material-ui/core/ListItemIcon'
import ListItemText from '@material-ui/core/ListItemText'
import BackupIcon from '@material-ui/icons/Backup'
import SearchIcon from '@material-ui/icons/Search'
import AccountIcon from '@material-ui/icons/AccountCircle'
import DocumentationIcon from '@material-ui/icons/Help'
import HomeIcon from '@material-ui/icons/Home'
import ArchiveIcon from '@material-ui/icons/Storage'
import EncIcon from '@material-ui/icons/Assessment'
import DevelIcon from '@material-ui/icons/ReportProblem'
import { Link, withRouter } from 'react-router-dom'
import { compose } from 'recompose'
import { Avatar, MuiThemeProvider } from '@material-ui/core'
import { genTheme, repoTheme, archiveTheme, encTheme, appBase } from '../config'
import { ErrorSnacks } from './errors'

const drawerWidth = 200

const toolbarTitles = {
  '/': 'Welcome',
  '/repo': 'Raw Code Outputs',
  '/upload': 'Upload Your Own Data',
  '/profile': 'Your Profile',
  '/docs': 'Documentation',
  '/archive': 'Code Independent Data',
  '/enc': 'The Material Perspective',
  '/dev': 'Developer and Operator Functions'
}

const toolbarThemes = {
  '/': genTheme,
  '/repo': repoTheme,
  '/upload': repoTheme,
  '/profile': genTheme,
  '/docs': genTheme,
  '/archive': archiveTheme,
  '/enc': encTheme,
  '/dev': genTheme
}

class Navigation extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    children: PropTypes.any,
    location: PropTypes.object.isRequired
  }

  static styles = theme => ({
    root: {
      flexGrow: 1
    },
    flex: {
      flexGrow: 1
    },
    appFrame: {
      zIndex: 1,
      overflow: 'hidden',
      position: 'relative',
      display: 'flex',
      width: '100%',
      height: '100vh'
    },
    appBar: {
      width: `calc(100% - ${drawerWidth}px)`,
      marginLeft: drawerWidth
    },
    drawerPaper: {
      position: 'relative',
      width: drawerWidth
    },
    content: {
      flexGrow: 1,
      backgroundColor: theme.palette.background.default,
      padding: theme.spacing.unit * 3,
      width: '100%',
      overflow: 'scroll'
    },
    toolbar: theme.mixins.toolbar,
    link: {
      textDecoration: 'none',
      color: theme.palette.text.primary
    }
  });

  renderTitle() {
    return (
      <Toolbar>
        <Typography variant="title" color="inherit" noWrap>
          nomad xt
        </Typography>
      </Toolbar>
    )
  }

  renderMenu() {
    const {pathname} = this.props.location

    return (
      <div>
        <MenuList>
          <MenuItem component={Link} to="/" selected={ pathname === '/' }>
            <ListItemIcon>
              <HomeIcon />
            </ListItemIcon>
            <ListItemText inset primary="Home"/>
          </MenuItem>
          <MenuItem component={Link} to="/repo" selected={ pathname.startsWith('/repo') }>
            <ListItemIcon>
              <SearchIcon style={{fill: repoTheme.palette.primary.main}}/>
            </ListItemIcon>
            <ListItemText inset primary="Repository"/>
          </MenuItem>
          <MenuItem component={Link} to="/upload" selected={ pathname === '/upload' }>
            <ListItemIcon>
              <BackupIcon style={{fill: repoTheme.palette.primary.main}}/>
            </ListItemIcon>
            <ListItemText inset primary="Upload"/>
          </MenuItem>
          <MenuItem component={Link} to="/archive" selected={ pathname.startsWith('/archive') }>
            <ListItemIcon>
              <ArchiveIcon style={{fill: archiveTheme.palette.primary.main}}/>
            </ListItemIcon>
            <ListItemText inset primary="Archive"/>
          </MenuItem>
          <MenuItem component={Link} to="/enc" selected={ pathname.startsWith('/enc') }>
            <ListItemIcon>
              <EncIcon style={{fill: encTheme.palette.primary.main}}/>
            </ListItemIcon>
            <ListItemText inset primary="Encyclopedia"/>
          </MenuItem>
        </MenuList>
        <Divider/>
        <MenuList>
          <MenuItem component={Link} to="/profile" selected={ pathname === '/profile' }>
            <ListItemIcon>
              <AccountIcon />
            </ListItemIcon>
            <ListItemText inset primary="Profil"/>
          </MenuItem>
          <MenuItem component={Link} to="/docs" selected={ pathname === '/docs' }>
            <ListItemIcon>
              <DocumentationIcon />
            </ListItemIcon>
            <ListItemText inset primary="Documentation"/>
          </MenuItem>
          <MenuItem component={Link} to="/dev" selected={ pathname === '/dev' }>
            <ListItemIcon>
              <DevelIcon />
            </ListItemIcon>
            <ListItemText inset primary="Development"/>
          </MenuItem>
        </MenuList>
      </div>
    )
  }

  render() {
    const { classes, children, location: { pathname } } = this.props

    const drawer = (
      <Drawer variant="permanent" classes={{ paper: classes.drawerPaper }} anchor="left">
        <div className={classes.toolbar}>
          {this.renderTitle()}
        </div>
        <Divider />
        {this.renderMenu()}
      </Drawer>
    )

    const selected = dct => {
      const key = Object.keys(dct).find(key => {
        return key === pathname || (key.length > 1 && pathname.startsWith(key))
      })
      return dct[key]
    }
    const theme = selected(toolbarThemes)

    return (
      <div className={classes.root}>
        <div className={classes.appFrame}>
          <MuiThemeProvider theme={theme}>
            <AppBar
              position="absolute"
              className={classes.appBar}
            >
              <Toolbar>
                <Typography variant="title" color="inherit" noWrap className={classes.flex}>
                  {selected(toolbarTitles)}
                </Typography>
                <Avatar src={`${appBase}/me.jpg`}/>
              </Toolbar>
            </AppBar>
          </MuiThemeProvider>
          {drawer}
          <MuiThemeProvider theme={theme}>
            <main className={classes.content}>
              <div className={classes.toolbar} />
              <ErrorSnacks>
                {children}
              </ErrorSnacks>
            </main>
          </MuiThemeProvider>
        </div>
      </div>
    )
  }
}

export default compose(withRouter, withStyles(Navigation.styles))(Navigation)
