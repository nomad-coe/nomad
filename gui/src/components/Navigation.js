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
import AnalyticsIcon from '@material-ui/icons/Settings'
import DevelIcon from '@material-ui/icons/ReportProblem'
import ChevronLeftIcon from '@material-ui/icons/ChevronLeft'
import MenuIcon from '@material-ui/icons/Menu'
import { Link, withRouter } from 'react-router-dom'
import { compose } from 'recompose'
import { MuiThemeProvider, IconButton, Button, Checkbox, FormLabel, DialogTitle, DialogContent, DialogContentText, TextField, DialogActions, Dialog } from '@material-ui/core'
import { genTheme, repoTheme, archiveTheme, encTheme, analyticsTheme } from '../config'
import { ErrorSnacks } from './errors'
import classNames from 'classnames'
import { HelpContext } from './help'
import { withApi } from './api'

const drawerWidth = 200

const toolbarTitles = {
  '/': 'Welcome',
  '/repo': 'Raw Code Outputs',
  '/upload': 'Upload Your Own Data',
  '/profile': 'Your Profile',
  '/docs': 'Documentation',
  '/archive': 'Code Independent Data',
  '/analytics': 'Big Data Analytics',
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
  '/analytics': analyticsTheme,
  '/dev': genTheme
}

class LoginLogoutComponent extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    userName: PropTypes.string,
    login: PropTypes.func.isRequired,
    logout: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {
      display: 'flex',
      alignItems: 'center',
      '& p': {
        marginRight: theme.spacing.unit * 2
      },
      '& button': {
        borderColor: theme.palette.getContrastText(theme.palette.primary.main),
        marginRight: theme.spacing.unit * 4
      }
    }
  })

  constructor(props) {
    super(props)
    this.handleLoginDialogClosed = this.handleLoginDialogClosed.bind(this)
    this.handleLogout = this.handleLogout.bind(this)
    this.handleChange = this.handleChange.bind(this)
  }

  state = {
    loginDialogOpen: false,
    userName: null,
    password: null
  }

  handleLoginDialogClosed(withLogin) {
    this.setState({loginDialogOpen: false})
    if (withLogin) {
      this.props.login(this.state.userName, this.state.password)
    }
  }

  handleChange = name => event => {
    this.setState({
      [name]: event.target.value
    })
  }

  handleLogout() {
    this.props.logout()
  }

  render() {
    const { classes, userName } = this.props
    if (userName) {
      return (
        <div className={classes.root}>
          <Typography color="inherit" variant="body1">
            Welcome, {userName}
          </Typography>
          <Button color="inherit" variant="outlined" onClick={this.handleLogout}>Logout</Button>
        </div>
      )
    } else {
      return (
        <div className={classes.root}>
          <Button color="inherit" variant="outlined" onClick={() => this.setState({loginDialogOpen: true})}>Login</Button>
          <Dialog
            open={this.state.loginDialogOpen}
            onClose={this.handleLoginDialogClosed}
          >
            <DialogTitle>Login</DialogTitle>
            <DialogContent>
              <DialogContentText>
                To login, please enter your email address and password. If you
                do not have an account, please go to the nomad repository and
                create one.
              </DialogContentText>
              <TextField
                autoFocus
                margin="dense"
                id="uaseName"
                label="Email Address"
                type="email"
                fullWidth
                value={this.state.userName}
                onChange={this.handleChange('userName')}
              />
              <TextField
                margin="dense"
                id="password"
                label="Password"
                type="password"
                fullWidth
                value={this.state.password}
                onChange={this.handleChange('password')}
              />
            </DialogContent>
            <DialogActions>
              <Button onClick={this.handleLoginDialogClosed} color="primary">
                Cancel
              </Button>
              <Button onClick={() => this.handleLoginDialogClosed(true)} color="primary">
                Login
              </Button>
            </DialogActions>
          </Dialog>
        </div>
      )
    }
  }
}

const LoginLogout = compose(withApi(false), withStyles(LoginLogoutComponent.styles))(LoginLogoutComponent)

class Navigation extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    children: PropTypes.any,
    location: PropTypes.object.isRequired
  }

  static styles = theme => ({
    // root: {
    //   flexGrow: 1
    // },
    flex: {
      flexGrow: 1
    },
    root: {
      flexGrow: 1
      // height: 440,
      // zIndex: 1,
      // overflow: 'hidden',
      // position: 'relative',
      // display: 'flex'
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
      zIndex: theme.zIndex.drawer + 1,
      paddingRight: theme.spacing.unit * 3,
      transition: theme.transitions.create(['width', 'margin'], {
        easing: theme.transitions.easing.sharp,
        duration: theme.transitions.duration.leavingScreen
      })
    },
    appBarShift: {
      marginLeft: drawerWidth,
      paddingRight: theme.spacing.unit * 0,
      width: `calc(100% - ${drawerWidth}px)`,
      transition: theme.transitions.create(['width', 'margin'], {
        easing: theme.transitions.easing.sharp,
        duration: theme.transitions.duration.enteringScreen
      })
    },
    menuButton: {
      marginLeft: 12,
      marginRight: 36
    },
    hide: {
      display: 'none'
    },
    drawerPaper: {
      position: 'relative',
      whiteSpace: 'nowrap',
      width: drawerWidth,
      transition: theme.transitions.create('width', {
        easing: theme.transitions.easing.sharp,
        duration: theme.transitions.duration.enteringScreen
      })
    },
    drawerPaperClose: {
      overflowX: 'hidden',
      transition: theme.transitions.create('width', {
        easing: theme.transitions.easing.sharp,
        duration: theme.transitions.duration.leavingScreen
      }),
      width: theme.spacing.unit * 7,
      [theme.breakpoints.up('sm')]: {
        width: theme.spacing.unit * 9
      }
    },
    toolbar: {
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'flex-end',
      padding: '0 8px',
      ...theme.mixins.toolbar
    },
    // content: {
    //   flexGrow: 1,
    //   backgroundColor: theme.palette.background.default,
    //   padding: theme.spacing.unit * 3,
    // },
    content: {
      flexGrow: 1,
      backgroundColor: theme.palette.background.default,
      padding: theme.spacing.unit * 3,
      width: '100%',
      overflow: 'scroll'
    },
    link: {
      textDecoration: 'none',
      color: theme.palette.text.primary
    },
    menuItem: {
      paddingLeft: theme.spacing.unit * 3
    },
    barActions: {
      display: 'flex',
      alignItems: 'center'
    },
    barSelect: {
      color: `${theme.palette.getContrastText(theme.palette.primary.main)} !important`
    }
  })

  state = {
    open: false
  };

  handleDrawerOpen = () => {
    this.setState({ open: true })
  }

  handleDrawerClose = () => {
    this.setState({ open: false })
  }

  renderTitle() {
    return (
      <Toolbar>
        <Typography style={{fontSize: 24}} color="inherit" noWrap>
          nomad <strong>xt</strong>
        </Typography>
      </Toolbar>
    )
  }

  renderMenu() {
    const {classes} = this.props
    const {pathname} = this.props.location

    return (
      <div>
        <MenuList>
          <MenuItem className={classes.menuItem} component={Link} to="/" selected={ pathname === '/' }>
            <ListItemIcon>
              <HomeIcon />
            </ListItemIcon>
            <ListItemText inset primary="Home"/>
          </MenuItem>
          <MenuItem className={classes.menuItem} component={Link} to="/repo" selected={ pathname.startsWith('/repo') }>
            <ListItemIcon>
              <SearchIcon style={{fill: repoTheme.palette.primary.main}}/>
            </ListItemIcon>
            <ListItemText inset primary="Repository"/>
          </MenuItem>
          <MenuItem className={classes.menuItem} component={Link} to="/upload" selected={ pathname === '/upload' }>
            <ListItemIcon>
              <BackupIcon style={{fill: repoTheme.palette.primary.main}}/>
            </ListItemIcon>
            <ListItemText inset primary="Upload"/>
          </MenuItem>
          <MenuItem className={classes.menuItem} component={Link} to="/archive" selected={ pathname.startsWith('/archive') }>
            <ListItemIcon>
              <ArchiveIcon style={{fill: archiveTheme.palette.primary.main}}/>
            </ListItemIcon>
            <ListItemText inset primary="Archive"/>
          </MenuItem>
          <MenuItem className={classes.menuItem} component={Link} to="/enc" selected={ pathname.startsWith('/enc') }>
            <ListItemIcon>
              <EncIcon style={{fill: encTheme.palette.primary.main}}/>
            </ListItemIcon>
            <ListItemText inset primary="Encyclopedia"/>
          </MenuItem>
          <MenuItem className={classes.menuItem} component={Link} to="/analytics" selected={ pathname.startsWith('/analytics') }>
            <ListItemIcon>
              <AnalyticsIcon style={{fill: analyticsTheme.palette.primary.main}}/>
            </ListItemIcon>
            <ListItemText inset primary="Analytics"/>
          </MenuItem>
        </MenuList>
        <Divider/>
        <MenuList>
          <MenuItem className={classes.menuItem} component={Link} to="/profile" selected={ pathname === '/profile' }>
            <ListItemIcon>
              <AccountIcon />
            </ListItemIcon>
            <ListItemText inset primary="Profil"/>
          </MenuItem>
          <MenuItem className={classes.menuItem} component={Link} to="/docs" selected={ pathname === '/docs' }>
            <ListItemIcon>
              <DocumentationIcon />
            </ListItemIcon>
            <ListItemText inset primary="Documentation"/>
          </MenuItem>
          <MenuItem className={classes.menuItem} component={Link} to="/dev" selected={ pathname === '/dev' }>
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
      <Drawer variant="permanent"
        open={this.state.open}
        classes={{ paper: classNames(classes.drawerPaper, !this.state.open && classes.drawerPaperClose) }}
        anchor="left"
      >
        <div className={classes.toolbar}>
          <IconButton onClick={this.handleDrawerClose}>
            <ChevronLeftIcon/>
          </IconButton>
        </div>
        <Divider />
        {/* <div className={classes.toolbar}>
          {this.renderTitle()}
        </div> */}
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
              className={classNames(classes.appBar, this.state.open && classes.appBarShift)}
            >
              <Toolbar disableGutters={!this.state.open}>
                <IconButton
                  color="inherit"
                  aria-label="Open drawer"
                  onClick={this.handleDrawerOpen}
                  className={classNames(classes.menuButton, this.state.open && classes.hide)}
                >
                  <MenuIcon />
                </IconButton>
                <Typography variant="h6" color="inherit" noWrap className={classes.flex}>
                  {selected(toolbarTitles)}
                </Typography>
                <div className={classes.barActions}>
                  <LoginLogout />
                  <FormLabel className={classes.barSelect} >Show help</FormLabel>
                  <HelpContext.Consumer>{
                    help => (
                      <Checkbox
                        checked={!help.someClosed()} indeterminate={!help.allClosed() && help.someClosed()}
                        onClick={() => help.switchHelp()}
                        classes={{root: classes.barSelect, checked: classes.barSelect}} />
                    )
                  }</HelpContext.Consumer>
                </div>
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
