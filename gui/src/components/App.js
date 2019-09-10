// trigger rebuild

import React from 'react'
import PropTypes, { instanceOf } from 'prop-types'
import { compose } from 'recompose'
import classNames from 'classnames'
import { MuiThemeProvider, withStyles } from '@material-ui/core/styles'
import { IconButton, LinearProgress, ListItemIcon, ListItemText,
  MenuList, MenuItem, Typography, Drawer, AppBar, Toolbar, Divider, Button, DialogContent, DialogTitle, DialogActions, Dialog, Tooltip } from '@material-ui/core'
import { Switch, Route, Link, withRouter } from 'react-router-dom'
import BackupIcon from '@material-ui/icons/Backup'
import SearchIcon from '@material-ui/icons/Search'
import AboutIcon from '@material-ui/icons/Help'
import MetainfoIcon from '@material-ui/icons/Info'
import ChevronLeftIcon from '@material-ui/icons/ChevronLeft'
import MenuIcon from '@material-ui/icons/Menu'
import {help as searchHelp, default as SearchPage} from './search/SearchPage'
import HelpDialog from './Help'
import { ApiProvider, withApi } from './api'
import { ErrorSnacks, withErrors } from './errors'
import Calc from './entry/Calc'
import About from './About'
import LoginLogout from './LoginLogout'
import { genTheme, repoTheme, archiveTheme, appBase } from '../config'
import { DomainProvider, withDomain } from './domains'
import {help as metainfoHelp, default as MetaInfoBrowser} from './metaInfoBrowser/MetaInfoBrowser'
import packageJson from '../../package.json'
import { Cookies, withCookies } from 'react-cookie'
import Markdown from './Markdown'
import {help as uploadHelp, default as Uploads} from './uploads/Uploads'
import ResolvePID from './entry/ResolvePID'
import DatasetPage from './DatasetPage'
import { capitalize } from '../utils'

export class VersionMismatch extends Error {
  constructor(msg) {
    super(msg)
    this.name = 'VersionMismatch'
  }
}

const drawerWidth = 200



class NavigationUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    children: PropTypes.any,
    location: PropTypes.object.isRequired,
    loading: PropTypes.number.isRequired,
    raiseError: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {},
    title: {
      marginLeft: theme.spacing.unit,
      flexGrow: 1,
      display: 'flex',
      alignItems: 'center',
      alignContent: 'flex-start'
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
      marginLeft: theme.spacing.unit
    },
    helpButton: {
      marginLeft: theme.spacing.unit
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
      marginRight: 0
    },
    barButtonDisabled: {
      marginRight: 0
    }
  })

  constructor(props) {
    super(props)

    this.state = {
      open: false
    }
  }

  toolbarTitles = {
    '/': 'About, Documentation, Getting Help',
    '/search': 'Find and Download Data',
    '/uploads': 'Upload and Publish Data',
    '/metainfo': 'The NOMAD Meta Info',
    '/entry': capitalize(this.props.domain.entryLabel),
    '/dataset': 'Dataset'
  }

  toolbarThemes = {
    '/': genTheme,
    '/search': repoTheme,
    '/uploads': repoTheme,
    '/entry': repoTheme,
    '/dataset': repoTheme,
    '/metainfo': archiveTheme
  }

  toolbarHelp = {
    '/': null,
    '/search': {title: 'How to find and download data', content: searchHelp},
    '/uploads': {title: 'How to upload data', content: uploadHelp},
    '/metainfo': {title: 'About the NOMAD meta-info', content: metainfoHelp}
  }

  componentDidMount() {
    fetch(`${appBase}/meta.json`)
      .then((response) => response.json())
      .then((meta) => {
        if (meta.version !== packageJson.version) {
          // this should not happen, if we setup the web servers correctly
          console.error('GUI API version mismatch')
        }
      })
      .catch(() => {
        console.log('Could not validate version, continue...')
      })
  }

  handleDrawerEvent(isOpen) {
    this.setState({ open: !isOpen, openIsSet: true })
  }

  render() {
    const { classes, children, location: { pathname }, loading } = this.props
    const { toolbarThemes, toolbarHelp, toolbarTitles } = this

    const selected = dct => {
      const key = Object.keys(dct).find(key => {
        return key === pathname || (key.length > 1 && pathname.startsWith(key))
      })
      return dct[key]
    }

    const theme = selected(toolbarThemes)
    const help = selected(toolbarHelp)

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
                  onClick={() => this.handleDrawerEvent(this.state.open)}
                  className={classNames(classes.menuButton, this.state.open && classes.hide)}
                >
                  <MenuIcon />
                </IconButton>
                <div className={classes.title}>
                  <Typography variant="h6" color="inherit" noWrap>
                    {selected(toolbarTitles)}
                  </Typography>
                  {help ? <HelpDialog color="inherit" maxWidth="md" classes={{root: classes.helpButton}} {...help}/> : ''}
                </div>
                <div className={classes.barActions}>
                  <LoginLogout variant="outlined" color="inherit" classes={{button: classes.barButton, buttonDisabled: classes.barButtonDisabled}} />
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
                <IconButton onClick={() => this.handleDrawerEvent(this.state.open)}>
                  <ChevronLeftIcon/>
                </IconButton>
              </div>

              <MenuList>
                <Tooltip title="Upload and publish data">
                  <MenuItem className={classes.menuItem} component={Link} to="/uploads" selected={ pathname === '/uploads' }>
                    <ListItemIcon>
                      <BackupIcon style={{fill: repoTheme.palette.primary.main}}/>
                    </ListItemIcon>
                    <ListItemText inset primary="Upload"/>
                  </MenuItem>
                </Tooltip>
                <Tooltip title="Find and download data">
                  <MenuItem className={classes.menuItem} component={Link} to="/search" selected={ pathname.startsWith('/repo') }>
                    <ListItemIcon>
                      <SearchIcon style={{fill: repoTheme.palette.primary.main}}/>
                    </ListItemIcon>
                    <ListItemText inset primary="Search"/>
                  </MenuItem>
                </Tooltip>
                <Divider />
                <Tooltip title="Browse the archive schema">
                  <MenuItem className={classes.menuItem} component={Link} to="/metainfo" selected={ pathname === '/metainfo' }>
                    <ListItemIcon>
                      <MetainfoIcon style={{fill: archiveTheme.palette.primary.main}}/>
                    </ListItemIcon>
                    <ListItemText inset primary="Meta Info"/>
                  </MenuItem>
                </Tooltip>
                <Divider />
                <Tooltip title="About, Documentation, Getting Help">
                  <MenuItem className={classes.menuItem} component={Link} to="/" selected={ pathname === '/' }>
                    <ListItemIcon>
                      <AboutIcon />
                    </ListItemIcon>
                    <ListItemText inset primary="Help"/>
                  </MenuItem>
                </Tooltip>
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

const Navigation = compose(withRouter, withErrors, withApi(false), withDomain, withStyles(NavigationUnstyled.styles))(NavigationUnstyled)

class LicenseAgreementUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    cookies: instanceOf(Cookies).isRequired
  }

  static styles = theme => ({
    content: {
      backgroundColor: theme.palette.primary.main
    },
    button: {
      color: 'white'
    }
  })

  constructor(props) {
    super(props)

    this.handleClosed = this.handleClosed.bind(this)
  }

  state = {
    accepted: this.props.cookies.get('terms-accepted')
  }

  handleClosed(accepted) {
    if (accepted) {
      this.props.cookies.set('terms-accepted', true)
      this.setState({accepted: true})
    }
  }

  render() {
    return (
      <div>
        <Dialog
          disableBackdropClick disableEscapeKeyDown
          open={!this.state.accepted}
        >
          <DialogTitle>Terms of Use</DialogTitle>
          <DialogContent>
            <Markdown>{`
              By uploading and downloading data, you agree to the
              [terms of use](https://www.nomad-coe.eu/the-project/nomad-repository/nomad-repository-terms).

              Uploaded data is licensed under the Creative Commons Attribution license
              ([CC BY 3.0](https://creativecommons.org/licenses/by/3.0/)). You can publish
              uploaded data with an *embargo*. Data with an *embargo* is only visible to
              you and users you share your data with. The *embargo period* lasts up to 36 month.
              After the *embargo* your published data will be public. **Note that public data
              is visible to others and files become downloadable by everyone.**

              This web-site uses *cookies*. By using this web-site you agree to our use
              of *cookies*. [Learn more](https://www.cookiesandyou.com/).
              `}
            </Markdown>
          </DialogContent>
          <DialogActions>
            <Button onClick={() => this.handleClosed(true)} color="primary">
              Accept
            </Button>
          </DialogActions>
        </Dialog>
      </div>
    )
  }
}

const LicenseAgreement = compose(withCookies, withStyles(LicenseAgreementUnstyled.styles))(LicenseAgreementUnstyled)

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
    'entry': {
      path: '/entry/id/:uploadId/:calcId',
      key: (props) => `entry/id/${props.match.params.uploadId}/${props.match.params.uploadId}`,
      render: props => {
        const { match, ...rest } = props
        if (match && match.params.uploadId && match.params.calcId) {
          return (<Calc {...rest} uploadId={match.params.uploadId} calcId={match.params.calcId} />)
        } else {
          return ''
        }
      }
    },
    'dataset': {
      path: '/dataset/id/:datasetId',
      key: (props) => `dataset/id/${props.match.params.datasetId}`,
      render: props => {
        const { match, ...rest } = props
        if (match && match.params.datasetId) {
          return (<DatasetPage {...rest} datasetId={match.params.datasetId} />)
        } else {
          return ''
        }
      }
    },
    'entry_pid': {
      path: '/entry/pid/:pid',
      key: (props) => `entry/pid/${props.match.params.pid}`,
      render: props => {
        const { match, ...rest } = props
        if (match && match.params.pid) {
          return (<ResolvePID {...rest} pid={match.params.pid} />)
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
        </ErrorSnacks>
        <LicenseAgreement />
      </MuiThemeProvider>
    )
  }
}
