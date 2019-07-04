// trigger rebuild

import React from 'react'
import PropTypes from 'prop-types'
import { compose } from 'recompose'
import classNames from 'classnames'
import { MuiThemeProvider, withStyles } from '@material-ui/core/styles'
import { IconButton, Checkbox, FormLabel, LinearProgress, ListItemIcon, ListItemText,
  MenuList, MenuItem, Typography, Drawer, AppBar, Toolbar, Divider } from '@material-ui/core'
import { BrowserRouter, Switch, Route, Link, withRouter } from 'react-router-dom'
import BackupIcon from '@material-ui/icons/Backup'
import SearchIcon from '@material-ui/icons/Search'
import AboutIcon from '@material-ui/icons/Home'
import MetainfoIcon from '@material-ui/icons/Info'
import ChevronLeftIcon from '@material-ui/icons/ChevronLeft'
import MenuIcon from '@material-ui/icons/Menu'

import Uploads from './uploads/Uploads'
import SearchPage from './search/SearchPage'
import { HelpProvider, HelpContext } from './help'
import { ApiProvider, withApi } from './api'
import { ErrorSnacks } from './errors'
import Calc from './entry/Calc'
import About from './About'
import LoginLogout from './LoginLogout'
import { genTheme, repoTheme, archiveTheme } from '../config'
import { DomainProvider } from './domains'
import MetaInfoBrowser from './metaInfoBrowser/MetaInfoBrowser'
import packageJson from '../../package.json'

const drawerWidth = 200

const toolbarTitles = {
  '/': 'About nomad@FAIRDI',
  '/search': 'Search',
  '/uploads': 'Upload Your Own Data',
  '/metainfo': 'The Nomad Meta Info'
}

const toolbarThemes = {
  '/': genTheme,
  '/search': repoTheme,
  '/uploads': repoTheme,
  '/metainfo': archiveTheme
}

class NavigationUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    children: PropTypes.any,
    location: PropTypes.object.isRequired,
    loading: PropTypes.number.isRequired
  }

  static styles = theme => ({
    root: {},
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
    content: {
      flexGrow: 1,
      backgroundColor: theme.palette.background.default,
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
    },
    barButton: {
      borderColor: theme.palette.getContrastText(theme.palette.primary.main),
      marginRight: theme.spacing.unit * 4
    },
    barButtonDisabled: {
      marginRight: theme.spacing.unit * 4
    }
  })

  state = {
    open: false
  }

  constructor(props) {
    super(props)

    this.handleDrawerOpen = this.handleDrawerOpen.bind(this)
    this.handleDrawerClose = this.handleDrawerClose.bind(this)
  }

  componentDidMount() {
    fetch('/meta.json')
      .then((response) => response.json())
      .then((meta) => {
        if (meta.version !== packageJson.version) {
          console.log('Different version, hard reloading...')
          window.location.reload(true)
        }
      })
      .catch(() => {
        console.log('Could not validate version, continue...')
      })
  }

  handleDrawerOpen() {
    this.setState({ open: true })
  }

  handleDrawerClose() {
    this.setState({ open: false })
  }

  render() {
    const { classes, children, location: { pathname }, loading } = this.props

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
                  onClick={this.handleDrawerOpen}
                  className={classNames(classes.menuButton, this.state.open && classes.hide)}
                >
                  <MenuIcon />
                </IconButton>
                <Typography variant="h6" color="inherit" noWrap className={classes.flex}>
                  {selected(toolbarTitles)}
                </Typography>
                <div className={classes.barActions}>
                  <LoginLogout variant="outlined" color="inherit" classes={{button: classes.barButton, buttonDisabled: classes.barButtonDisabled}} />
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
              {loading ? <LinearProgress color="primary" /> : ''}
            </AppBar>

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

              <MenuList>
                <MenuItem className={classes.menuItem} component={Link} to="/" selected={ pathname === '/' }>
                  <ListItemIcon>
                    <AboutIcon />
                  </ListItemIcon>
                  <ListItemText inset primary="About"/>
                </MenuItem>
                <Divider />
                <MenuItem className={classes.menuItem} component={Link} to="/search" selected={ pathname.startsWith('/repo') }>
                  <ListItemIcon>
                    <SearchIcon style={{fill: repoTheme.palette.primary.main}}/>
                  </ListItemIcon>
                  <ListItemText inset primary="Search"/>
                </MenuItem>
                <MenuItem className={classes.menuItem} component={Link} to="/uploads" selected={ pathname === '/uploads' }>
                  <ListItemIcon>
                    <BackupIcon style={{fill: repoTheme.palette.primary.main}}/>
                  </ListItemIcon>
                  <ListItemText inset primary="Upload"/>
                </MenuItem>
                <Divider />
                <MenuItem className={classes.menuItem} component={Link} to="/metainfo" selected={ pathname === '/metainfo' }>
                  <ListItemIcon>
                    <MetainfoIcon style={{fill: archiveTheme.palette.primary.main}}/>
                  </ListItemIcon>
                  <ListItemText inset primary="Meta Info"/>
                </MenuItem>
              </MenuList>
            </Drawer>

            <main className={classes.content}>
              <div className={classes.toolbar} />
              {children}
            </main>

          </MuiThemeProvider>
        </div>
      </div>
    )
  }
}

const Navigation = compose(withRouter, withApi(false), withStyles(NavigationUnstyled.styles))(NavigationUnstyled)

export default class App extends React.Component {
  constructor(props) {
    super(props)
    this.renderChildren.bind(this)
  }

  routes = {
    'about': {
      exact: true,
      singleton: true,
      path: '/',
      render: props => <About {...props} />
    },
    'search': {
      exact: true,
      singleton: true,
      path: '/search',
      render: props => <SearchPage {...props} />
    },
    'searchEntry': {
      path: '/search/:uploadId/:calcId',
      key: (props) => `searchEntry/${props.match.params.uploadId}/${props.match.params.uploadId}`,
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
      singleton: true,
      path: '/uploads',
      render: props => <Uploads {...props} />
    },
    'uploadedEntry': {
      path: '/uploads/:uploadId/:calcId',
      key: (props) => `uploadedEntry/${props.match.params.uploadId}/${props.match.params.uploadId}`,
      render: props => {
        const { match, ...rest } = props
        if (match && match.params.uploadId && match.params.calcId) {
          return (<Calc {...rest} uploadId={match.params.uploadId} calcId={match.params.calcId} />)
        } else {
          return ''
        }
      }
    },
    'metainfo': {
      exact: true,
      path: '/metainfo',
      render: props => <MetaInfoBrowser {...props} />
    },
    'metainfoEntry': {
      path: '/metainfo/:metainfo',
      key: props => `metainfo/${props.match.params.metainfo}`,
      render: props => <MetaInfoBrowser metainfo={props.match.params.metainfo} {...props} />
    }
  }

  renderChildren(routeKey, props) {
    // const { match, ...rest } = props

    return (
      <div>
        {Object.keys(this.routes)
          .filter(route => this.routes[route].singleton || route === routeKey)
          .map(route => (
            <div
              key={route.key ? route.key(props) : route}
              style={{display: routeKey === route ? 'block' : 'none'}}
            >
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
          <BrowserRouter basename={process.env.PUBLIC_URL}>
            <HelpProvider>
              <ApiProvider>
                <DomainProvider>
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
                </DomainProvider>
              </ApiProvider>
            </HelpProvider>
          </BrowserRouter>
        </ErrorSnacks>
      </MuiThemeProvider>
    )
  }
}
