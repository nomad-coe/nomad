// trigger rebuild

import React from 'react'
import PropTypes, { instanceOf } from 'prop-types'
import { compose } from 'recompose'
import classNames from 'classnames'
import { MuiThemeProvider, withStyles } from '@material-ui/core/styles'
import { LinearProgress, ListItemIcon, ListItemText, MenuList, MenuItem, Typography,
  AppBar, Toolbar, Button, DialogContent, DialogTitle, DialogActions, Dialog, Tooltip,
  Snackbar, SnackbarContent } from '@material-ui/core'
import { Switch, Route, Link, withRouter } from 'react-router-dom'
import BackupIcon from '@material-ui/icons/Backup'
import SearchIcon from '@material-ui/icons/Search'
import UserDataIcon from '@material-ui/icons/AccountCircle'
import AboutIcon from '@material-ui/icons/Home'
import MetainfoIcon from '@material-ui/icons/Info'
import {help as searchHelp, default as SearchPage} from './search/SearchPage'
import HelpDialog from './Help'
import { ApiProvider, withApi } from './api'
import { ErrorSnacks, withErrors } from './errors'
import { help as entryHelp, default as EntryPage } from './entry/EntryPage'
import About from './About'
import LoginLogout from './LoginLogout'
import { guiBase, consent, nomadTheme } from '../config'
import { DomainProvider, withDomain } from './domains'
import {help as metainfoHelp, default as MetaInfoBrowser} from './metaInfoBrowser/MetaInfoBrowser'
import packageJson from '../../package.json'
import { Cookies, withCookies } from 'react-cookie'
import Markdown from './Markdown'
import {help as uploadHelp, default as Uploads} from './uploads/Uploads'
import ResolvePID from './entry/ResolvePID'
import DatasetPage from './DatasetPage'
import { capitalize } from '../utils'
import { amber } from '@material-ui/core/colors'
import KeepState from './KeepState'
import {help as userdataHelp, default as UserdataPage} from './UserdataPage'
import ResolveDOI from './dataset/ResolveDOI'

export const ScrollContext = React.createContext({scrollParentRef: null});

export class VersionMismatch extends Error {
  constructor(msg) {
    super(msg)
    this.name = 'VersionMismatch'
  }
}

function ReloadSnack() {
  return <Snackbar
    anchorOrigin={{
      vertical: 'bottom',
      horizontal: 'left'
    }}
    open
  >
    <SnackbarContent
      style={{backgroundColor: amber[700]}}
      message={<span>There is a new NOMAD version. Please press your browser's reload (or even shift+reload) button.</span>}
    />
  </Snackbar>
}

class NavigationUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    children: PropTypes.any,
    location: PropTypes.object.isRequired,
    loading: PropTypes.number.isRequired,
    raiseError: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {
      minWidth: 1024,
    },
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
      backgroundColor: '#20335D'
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
    toolbar: {
      paddingRight: theme.spacing.unit * 3
    },
    logo: {
      height: theme.spacing.unit * 7,
      marginRight: theme.spacing.unit * 2
    },
    menu: {
      display: 'inline-flex',
      padding: 0,
      width: '100%',
      backgroundColor: 'white'
    },
    content: {
      marginTop: theme.spacing.unit * 13,
      flexGrow: 1,
      backgroundColor: theme.palette.background.default,
      width: '100%',
      overflow: 'auto'
    },
    link: {
      textDecoration: 'none',
      color: theme.palette.text.primary
    },
    menuItemIcon: {
      marginRight: 0
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
    divider: {
      flexGrow: 1
    }
  })

  constructor(props) {
    super(props)
    this.scroll = {
      scrollParentRef: null
    }
    this.state = {
      open: false
    }
  }

  toolbarTitles = {
    '/': 'About, Documentation, Getting Help',
    '/search': 'Find and Download Data',
    '/uploads': 'Upload and Publish Data',
    '/userdata': 'Manage Your Data',
    '/metainfo': 'The NOMAD Meta Info',
    '/entry': capitalize(this.props.domain.entryLabel),
    '/dataset': 'Dataset'
  }

  toolbarHelp = {
    '/': null,
    '/search': {title: 'How to find and download data', content: searchHelp},
    '/uploads': {title: 'How to upload data', content: uploadHelp},
    '/userdata': {title: 'How to manage your data', content: userdataHelp},
    '/metainfo': {title: 'About the NOMAD meta-info', content: metainfoHelp},
    '/entry': {title: 'The entry page', content: entryHelp}
  }

  componentDidMount() {
    fetch(`${guiBase}/meta.json`, {
      method: 'GET',
      cache: 'no-cache',
      headers: {
        'Pragma': 'no-cache',
        'Cache-Control': 'no-cache, no-store'
      }
    }).then((response) => response.json())
      .then((meta) => {
        if (meta.version !== packageJson.version) {
          console.log('GUI API version mismatch')
          this.setState({showReloadSnack: true})
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
    const { toolbarHelp, toolbarTitles } = this
    const { showReloadSnack } = this.state

    const selected = dct => {
      const key = Object.keys(dct).find(key => {
        return key === pathname || (key.length > 1 && pathname.startsWith(key))
      })
      return dct[key]
    }

    const theme = nomadTheme
    const help = selected(toolbarHelp)

    return (
      <div className={classes.root}>
        <div className={classes.appFrame}>
          <MuiThemeProvider theme={theme}>
            { showReloadSnack ? <ReloadSnack/> : ''}
            <AppBar
              // position="absolute"
              position="fixed"
              className={classNames(classes.appBar, this.state.open && classes.appBarShift)}
            >
              <Toolbar classes={{root: classes.toolbar}}
                disableGutters
                // disableGutters={!this.state.open}
              >
                {/* <IconButton
                  color="inherit"
                  onClick={() => this.handleDrawerEvent(this.state.open)}
                  className={classNames(classes.menuButton, this.state.open && classes.hide)}
                >
                  <MenuIcon />
                </IconButton> */}
                <div className={classes.title}>
                  <img alt="The NOMAD logo" className={classes.logo} src="/nomad.png"></img>
                  <Typography variant="h6" color="inherit" noWrap>
                    {selected(toolbarTitles)}
                  </Typography>
                  {help ? <HelpDialog color="inherit" maxWidth="md" classes={{root: classes.helpButton}} {...help}/> : ''}
                </div>
                <div className={classes.barActions}>
                  <LoginLogout variant="outlined" color="inherit" classes={{button: classes.barButton}} />
                </div>
              </Toolbar>
              <MenuList classes={{root: classes.menu}}>
                <Tooltip title="Find and download data">
                  <MenuItem component={Link} to="/search" selected={ pathname.startsWith('/search') } dense>
                    <ListItemIcon classes={{root: classes.menuItemIcon}}>
                      <SearchIcon />
                    </ListItemIcon>
                    <ListItemText inset primary="Search"/>
                  </MenuItem>
                </Tooltip>
                <Tooltip title="Upload and publish data">
                  <MenuItem component={Link} to="/uploads" selected={ pathname === '/uploads' } dense>
                    <ListItemIcon classes={{root: classes.menuItemIcon}}>
                      <BackupIcon />
                    </ListItemIcon>
                    <ListItemText inset primary="Upload"/>
                  </MenuItem>
                </Tooltip>
                <Tooltip title="Manage your data">
                  <MenuItem component={Link} to="/userdata" selected={ pathname.startsWith('/userdata') } dense>
                    <ListItemIcon classes={{root: classes.menuItemIcon}}>
                      <UserDataIcon />
                    </ListItemIcon>
                    <ListItemText inset primary="Your data"/>
                  </MenuItem>
                </Tooltip>
                <div className={classes.divider} />
                <Tooltip title="NOMAD Repository and Archive">
                  <MenuItem component={Link} to="/" selected={ pathname === '/' } dense>
                    <ListItemIcon classes={{root: classes.menuItemIcon}}>
                      <AboutIcon />
                    </ListItemIcon>
                    <ListItemText inset primary="Overview"/>
                  </MenuItem>
                </Tooltip>
                <Tooltip title="Browse the archive schema">
                  <MenuItem component={Link} to="/metainfo" selected={ pathname === '/metainfo' } dense>
                    <ListItemIcon classes={{root: classes.menuItemIcon}}>
                      <MetainfoIcon />
                    </ListItemIcon>
                    <ListItemText inset primary="Meta Info"/>
                  </MenuItem>
                </Tooltip>
              </MenuList>
              {loading ? <LinearProgress color="primary" /> : ''}
            </AppBar>

            <main className={classes.content} ref={(ref) => this.scroll.scrollParentRef = ref}>
              <ScrollContext.Provider value={this.scroll}>
                {children}
              </ScrollContext.Provider>
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
            <Markdown>{consent}</Markdown>
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
    'userdata': {
      exact: true,
      singleton: true,
      path: '/userdata',
      render: props => <UserdataPage {...props} />
    },
    'entry': {
      path: '/entry/id/:uploadId/:calcId',
      key: (props) => `entry/id/${props.match.params.uploadId}/${props.match.params.uploadId}`,
      render: props => {
        const { match, ...rest } = props
        if (match && match.params.uploadId && match.params.calcId) {
          return (<EntryPage {...rest} uploadId={match.params.uploadId} calcId={match.params.calcId} />)
        } else {
          return ''
        }
      }
    },
    'entry_query': {
      exact: true,
      path: '/entry/query',
      render: props => <EntryPage {...props} query />
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
    'dataset_doi': {
      path: '/dataset/doi/:doi*',
      key: (props) => `dataset/doi/${props.match.params.doi}`,
      render: props => {
        const { match, ...rest } = props
        if (match && match.params.doi) {
          return (<ResolveDOI {...rest} doi={match.params.doi} />)
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
      singleton: true,
      render: props => <MetaInfoBrowser {...props} />
    },
    'metainfoEntry': {
      path: '/metainfo/:metainfo',
      key: props => `metainfo/${props.match.params.metainfo}`,
      render: props => <MetaInfoBrowser metainfo={props.match.params.metainfo} {...props} />
    }
  }

  renderChildren(routeKey, props) {
    return (
      <React.Fragment>
        {Object.keys(this.routes).map(route => <KeepState key={route}
          visible={routeKey === route}
          render={(props) => this.routes[route].render(props)}
          {...props} />)}
      </React.Fragment>
    )
  }

  render() {
    return (
      <MuiThemeProvider theme={nomadTheme}>
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
