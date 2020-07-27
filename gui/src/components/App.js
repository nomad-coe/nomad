import React, { useEffect, useState, useContext, useCallback, useRef } from 'react'
import PropTypes from 'prop-types'
import { compose } from 'recompose'
import classNames from 'classnames'
import { MuiThemeProvider, withStyles, makeStyles } from '@material-ui/core/styles'
import { LinearProgress, MenuList, Typography,
  AppBar, Toolbar, Button, DialogContent, DialogTitle, DialogActions, Dialog, Tooltip,
  Snackbar, SnackbarContent, FormGroup, FormControlLabel, Switch, IconButton } from '@material-ui/core'
import { Route, Link, withRouter, useLocation } from 'react-router-dom'
import BackupIcon from '@material-ui/icons/Backup'
import SearchIcon from '@material-ui/icons/Search'
import UserDataIcon from '@material-ui/icons/AccountCircle'
import AboutIcon from '@material-ui/icons/Home'
import FAQIcon from '@material-ui/icons/QuestionAnswer'
import MetainfoIcon from '@material-ui/icons/Info'
import DocIcon from '@material-ui/icons/Help'
import CodeIcon from '@material-ui/icons/Code'
import TermsIcon from '@material-ui/icons/Assignment'
import UnderstoodIcon from '@material-ui/icons/Check'
import {help as searchHelp, default as SearchPage} from './search/SearchPage'
import HelpDialog from './Help'
import { ApiProvider, withApi, apiContext } from './api'
import { ErrorSnacks, withErrors } from './errors'
import { help as entryHelp, default as EntryPage } from './entry/EntryPage'
import About from './About'
import LoginLogout from './LoginLogout'
import { guiBase, consent, nomadTheme, appBase, version } from '../config'
import {help as metainfoHelp, default as MetaInfoBrowser} from './metaInfoBrowser/MetaInfoBrowser'
import packageJson from '../../package.json'
import {help as uploadHelp, default as UploadPage} from './uploads/UploadPage'
import ResolvePID from './entry/ResolvePID'
import DatasetPage from './DatasetPage'
import { amber } from '@material-ui/core/colors'
import {help as userdataHelp, default as UserdataPage} from './UserdataPage'
import ResolveDOI from './dataset/ResolveDOI'
import FAQ from './FAQ'
import EntryQuery from './entry/EntryQuery'
import {matomo} from '../index'
import { useCookies } from 'react-cookie'
import Markdown from './Markdown'

export const ScrollContext = React.createContext({scrollParentRef: null})

function LoadingIndicator() {
  const {api} = useContext(apiContext)
  const [loading, setLoading] = useState(0)
  const handleOnLoading = useCallback(loading => setLoading(loading), [setLoading])
  useEffect(() => {
    api.onLoading(handleOnLoading)
    return () => api.removeOnLoading(handleOnLoading)
  }, [api, handleOnLoading])

  return loading ? <LinearProgress color="secondary" /> : ''
}

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
      message={<span>There is a new NOMAD version. Please press your browser&apos;s reload (or even shift+reload) button.</span>}
    />
  </Snackbar>
}

const useMainMenuItemStyles = makeStyles(theme => ({
  button: {
    margin: theme.spacing(1)
  }
}))

function MainMenuItem({tooltip, title, path, href, onClick, icon}) {
  const {pathname} = useLocation()
  const classes = useMainMenuItemStyles()
  const selected = path === pathname || (path !== '/' && pathname.startsWith(path))
  const rest = path ? {to: path, component: Link} : {href: href}
  return <Tooltip title={tooltip}>
    <Button
      className={classes.button}
      color={selected ? 'primary' : 'default'}
      size="small"
      startIcon={icon}
      onClick={onClick}
      {...rest}
    >
      {title}
    </Button>
  </Tooltip>
}
MainMenuItem.propTypes = {
  'tooltip': PropTypes.string.isRequired,
  'title': PropTypes.string.isRequired,
  'path': PropTypes.string,
  'href': PropTypes.string,
  'onClick': PropTypes.func,
  'icon': PropTypes.element.isRequired
}

const useBetaSnackStyles = makeStyles(theme => ({
  root: {},
  snack: {
    backgroundColor: amber[700]
  }
}))
function BetaSnack() {
  const classes = useBetaSnackStyles()
  const [understood, setUnderstood] = useState(false)

  if (!version) {
    console.log.warning('no version data available')
    return ''
  }

  if (!version.isBeta) {
    return ''
  }

  return <Snackbar className={classes.root}
    anchorOrigin={{
      vertical: 'bottom',
      horizontal: 'left'
    }}
    open={!understood}
  >
    <SnackbarContent
      className={classes.snack}
      message={<span style={{color: 'white'}}>
        You are using a beta version of NOMAD ({version.label}). {
          version.usesBetaData ? 'This version is not using the official data. Everything you upload here, might get lost.' : ''
        } Click <a style={{color: 'white'}} href={version.officialUrl}>here for the official NOMAD version</a>.
      </span>}
      action={[
        <IconButton key={0} color="inherit" onClick={() => setUnderstood(true)}>
          <UnderstoodIcon />
        </IconButton>
      ]}
    />
  </Snackbar>
}

function Consent() {
  const [cookies, setCookie] = useCookies()
  const [accepted, setAccepted] = useState(cookies['terms-accepted'])
  const [optOut, setOptOut] = useState(cookies['tracking-enabled'] === 'false')

  useEffect(() => {
    if (!optOut) {
      matomo.push(['setConsentGiven'])
    } else {
      matomo.push(['requireConsent'])
    }
  })

  const handleClosed = accepted => {
    if (accepted) {
      setCookie('terms-accepted', true)
      setCookie('tracking-enabled', !optOut)
      setAccepted(true)
    }
  }
  const handleOpen = () => {
    setCookie('terms-accepted', false)
    setAccepted(false)
  }

  return (
    <React.Fragment>
      <MainMenuItem
        title="Terms"
        onClick={handleOpen}
        tooltip="NOMAD's terms"
        icon={<TermsIcon/>}
      />
      <Dialog
        disableBackdropClick disableEscapeKeyDown
        open={!accepted}
      >
        <DialogTitle>Terms of Use</DialogTitle>
        <DialogContent>
          <Markdown>{consent}</Markdown>
          <FormGroup>
            <FormControlLabel
              control={<Switch
                checked={optOut}
                onChange={(e) => {
                  setOptOut(!optOut)
                }}
                color="primary"
              />}
              label="Do not provide information about your use of NOMAD (opt-out)."
            />
          </FormGroup>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => handleClosed(true)} color="primary">
            Accept
          </Button>
        </DialogActions>
      </Dialog>
    </React.Fragment>
  )
}

const useMainMenuStyles = makeStyles(theme => ({
  root: {
    display: 'inline-flex',
    padding: 0,
    width: '100%',
    backgroundColor: 'white'
  },
  divider: {
    width: theme.spacing(3)
  }
}))

function MainMenu() {
  const classes = useMainMenuStyles()

  // We keep the URL of those path where components keep meaningful state in the URL.
  // If the menu is used to comeback, the old URL is used. Therefore, it appears as
  // if the same component instance with the same state is still there.
  const {pathname, search} = useLocation()
  const historyRef = useRef({
    search: '/search',
    userdata: '/userdata'
  })
  const history = {...historyRef.current}
  Object.keys(historyRef.current).forEach(key => {
    if (pathname.startsWith('/' + key)) {
      historyRef.current[key] = pathname + (search || '')
      history[key] = '/' + key
    }
  })

  return <MenuList classes={{root: classes.root}}>
    <MainMenuItem
      title="Search"
      path={history.search}
      tooltip="Find and download data"
      icon={<SearchIcon/>}
    />
    <MainMenuItem
      title="Upload"
      path="/uploads"
      tooltip="Upload and publish data"
      icon={<BackupIcon/>}
    />
    <MainMenuItem
      title="Your data"
      path={history.userdata}
      tooltip="Manage your data"
      icon={<UserDataIcon/>}
    />
    <MainMenuItem
      title="Meta Info"
      path="/metainfo"
      tooltip="Browse the archive schema"
      icon={<MetainfoIcon/>}
    />
    <div className={classes.divider} />
    <MainMenuItem
      title="About"
      path="/"
      tooltip="NOMAD Repository and Archive"
      icon={<AboutIcon/>}
    />
    <MainMenuItem
      title="FAQ"
      href="https://nomad-lab.eu/repository-archive-faqs"
      tooltip="Frequently Asked Questions (FAQ)"
      icon={<FAQIcon/>}
    />
    <MainMenuItem
      title="Docs"
      href={`${appBase}/docs/index.html`}
      tooltip="The NOMAD documentation"
      icon={<DocIcon/>}
    />
    <MainMenuItem
      title="Sources"
      href="https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR"
      tooltip="NOMAD's Gitlab project"
      icon={<CodeIcon/>}
    />
    <Consent />
  </MenuList>
}

class NavigationUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    children: PropTypes.any,
    location: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {
      minWidth: 1024
    },
    title: {
      marginLeft: theme.spacing(1),
      flexGrow: 1,
      display: 'flex',
      alignItems: 'center',
      alignContent: 'flex-start',
      color: theme.palette.primary.main
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
      backgroundColor: 'white'
    },
    menuButton: {
      marginLeft: theme.spacing(1)
    },
    helpButton: {
      marginLeft: theme.spacing(1)
    },
    hide: {
      display: 'none'
    },
    toolbar: {
      paddingRight: theme.spacing(3)
    },
    logo: {
      height: theme.spacing(7),
      marginRight: theme.spacing(2)
    },
    content: {
      marginTop: theme.spacing(13),
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
    '/faq': 'Frequently Asked Questions',
    '/search': 'Find and Download Data',
    '/uploads': 'Upload and Publish Data',
    '/userdata': 'Manage Your Data',
    '/metainfo': 'The NOMAD Meta Info',
    '/entry': 'Entry',
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
        if (meta.commit !== packageJson.commit) {
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
    const { classes, children, location: { pathname } } = this.props
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
              >
                <div className={classes.title}>
                  <a href="https://nomad-lab.eu">
                    <img alt="The NOMAD logo" className={classes.logo} src={`${guiBase}/nomad.png`}></img>
                  </a>
                  <Typography variant="h6" color="inherit" noWrap>
                    {selected(toolbarTitles)}
                  </Typography>
                  {help ? <HelpDialog color="inherit" maxWidth="md" classes={{root: classes.helpButton}} {...help}/> : ''}
                </div>
                <div className={classes.barActions}>
                  <LoginLogout color="primary" classes={{button: classes.barButton}} />
                </div>
              </Toolbar>
              <MainMenu />
              <LoadingIndicator />
            </AppBar>

            <main className={classes.content} ref={(ref) => { this.scroll.scrollParentRef = ref }}>
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

const Navigation = compose(withRouter, withErrors, withApi(false), withStyles(NavigationUnstyled.styles))(NavigationUnstyled)

const routes = {
  'about': {
    exact: true,
    path: '/',
    component: About
  },
  'faq': {
    exact: true,
    path: '/faq',
    component: FAQ
  },
  'search': {
    exact: true,
    path: '/search',
    component: SearchPage
  },
  'userdata': {
    exact: true,
    path: '/userdata',
    component: UserdataPage
  },
  'entry': {
    path: '/entry/id',
    component: EntryPage
  },
  'entry_query': {
    exact: true,
    path: '/entry/query',
    component: EntryQuery
  },
  'entry_pid': {
    path: '/entry/pid',
    component: ResolvePID
  },
  'dataset': {
    path: '/dataset/id',
    component: DatasetPage
  },
  'dataset_doi': {
    path: '/dataset/doi',
    component: ResolveDOI
  },
  'uploads': {
    exact: true,
    path: '/uploads',
    component: UploadPage
  },
  'metainfo': {
    path: '/metainfo',
    keepState: true,
    exists: false,
    component: MetaInfoBrowser
  }
}

class App extends React.PureComponent {
  render() {
    return (
      <MuiThemeProvider theme={nomadTheme}>
        <BetaSnack />
        <ErrorSnacks>
          <ApiProvider>
            <Navigation>
              {Object.keys(routes).map(routeKey => {
                const route = routes[routeKey]
                const { path, exact } = route
                return <Route key={routeKey} exact={exact} path={path}
                  // eslint-disable-next-line react/no-children-prop
                  children={props => {
                    if (route.keepState) {
                      if (props.match || route.exists) {
                        route.exists = true
                        return <route.component visible={props.match && true} {...props} />
                      } else {
                        return ''
                      }
                    } else {
                      return props.match && <route.component {...props} />
                    }
                  }}
                />
              })}
            </Navigation>
          </ApiProvider>
        </ErrorSnacks>
      </MuiThemeProvider>
    )
  }
}

const AppWithRouter = withRouter(App)
export default AppWithRouter
