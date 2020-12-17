/*
 * Copyright The NOMAD Authors.
 *
 * This file is part of NOMAD. See https://nomad-lab.eu for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
import React, { useEffect, useState, useContext, useCallback, useRef, useMemo } from 'react'
import PropTypes from 'prop-types'
import { compose } from 'recompose'
import classNames from 'classnames'
import { MuiThemeProvider, withStyles, makeStyles } from '@material-ui/core/styles'
import { LinearProgress, Typography,
  AppBar, Toolbar, Button, DialogContent, DialogTitle, DialogActions, Dialog,
  Snackbar, SnackbarContent, FormGroup, FormControlLabel, Switch, IconButton, Link as MuiLink } from '@material-ui/core'
import { Route, withRouter, useLocation } from 'react-router-dom'
import BackupIcon from '@material-ui/icons/Backup'
import SearchIcon from '@material-ui/icons/Search'
import UserDataIcon from '@material-ui/icons/AccountCircle'
import AboutIcon from '@material-ui/icons/Home'
import ForumIcon from '@material-ui/icons/QuestionAnswer'
import FAQIcon from '@material-ui/icons/LiveHelp'
import EncyclopediaIcon from '@material-ui/icons/Language'
import MetainfoIcon from '@material-ui/icons/Info'
import DocIcon from '@material-ui/icons/Help'
import CodeIcon from '@material-ui/icons/Code'
import TermsIcon from '@material-ui/icons/Assignment'
import UnderstoodIcon from '@material-ui/icons/Check'
import AnalyticsIcon from '@material-ui/icons/ShowChart'
import {help as searchHelp, default as SearchPage} from './search/SearchPage'
import HelpDialog from './Help'
import { ApiProvider, withApi, apiContext } from './api'
import { ErrorSnacks, withErrors } from './errors'
import { help as entryHelp, default as EntryPage } from './entry/EntryPage'
import About from './About'
import LoginLogout from './LoginLogout'
import { guiBase, consent, nomadTheme, appBase, version, oasis, aitoolkitEnabled, encyclopediaEnabled } from '../config'
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
import { help as metainfoHelp, MetainfoPage } from './archive/MetainfoBrowser'
import AIToolkitPage from './aitoolkit/AIToolkitPage'
import { MenuBarItem, MenuBar, MenuBarMenu } from './nav/MenuBar'

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
    console.warn('no version data available')
    return ''
  }

  if (!version.isBeta && !version.isTest) {
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
        You are using a {version.isBeta ? 'beta' : 'test'} version of NOMAD ({version.label}). {
          version.usesBetaData ? 'This version is not using the official data. Everything you upload here, might get lost.' : ''
        } Click <MuiLink style={{color: 'white'}} href={version.officialUrl}>here for the official NOMAD version</MuiLink>.
      </span>}
      action={[
        <IconButton key={0} color="inherit" onClick={() => setUnderstood(true)}>
          <UnderstoodIcon />
        </IconButton>
      ]}
    />
  </Snackbar>
}

function Consent(moreProps) {
  const [cookies, setCookie] = useCookies()
  const [accepted, setAccepted] = useState(cookies['terms-accepted'])
  const [optOut, setOptOut] = useState(cookies['tracking-enabled'] === 'false')
  const cookieOptions = useMemo(() => ({
    expires: new Date(2147483647 * 1000),
    path: '/' + guiBase.split('/').slice(1).join('/')
  }), [])

  useEffect(() => {
    if (!optOut) {
      matomo.push(['setConsentGiven'])
    } else {
      matomo.push(['requireConsent'])
    }
  })

  // Write again to push forwards Safari's hard-coded 7 days ITP window
  useEffect(() => {
    setCookie('terms-accepted', cookies['terms-accepted'], cookieOptions)
    setCookie('tracking-enabled', cookies['tracking-enabled'], cookieOptions)
  },
  // eslint-disable-next-line react-hooks/exhaustive-deps
  [])

  const handleClosed = accepted => {
    if (accepted) {
      setCookie('terms-accepted', true, cookieOptions)
      setCookie('tracking-enabled', !optOut, cookieOptions)
      setAccepted(true)
    }
  }
  const handleOpen = () => {
    setCookie('terms-accepted', false, cookieOptions)
    setAccepted(false)
  }

  return (
    <React.Fragment>
      <MenuBarItem
        name="terms"
        onClick={handleOpen}
        tooltip="The terms of service and cookie consent"
        icon={<TermsIcon/>}
        {...moreProps}
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

function MainMenu() {
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
  const route = Object.keys(routes).find(routeKey => pathname.startsWith(routes[routeKey].path))
  const routeNavPath = route && routes[route].navPath
  if (routeNavPath) {
    historyRef.current.navPath = routeNavPath
  }
  const selected = (route && routes[route].navPath) || historyRef.current.navPath || (route && routes[route].defaultNavPath) || 'publish/uploads'

  return <MenuBar selected={selected}>
    <MenuBarMenu name="publish" label="Publish" route="/uploads" icon={<BackupIcon/>}>
      <MenuBarItem
        name="uploads" label="Upload" route="/uploads" isDefault
        tooltip="Upload and publish new data" icon={<SearchIcon />}
      />
      <MenuBarItem
        label="Your data" name="userdata" route={history.userdata}
        tooltip="Manage your uploaded data" icon={<UserDataIcon />}
      />
    </MenuBarMenu>
    <MenuBarMenu name="explore" route={history.search} icon={<SearchIcon/>}>
      <MenuBarItem
        name="search" route={history.search}
        tooltip="Find and download data"
      />
      {encyclopediaEnabled && <MenuBarItem
        name="encyclopedia"
        href={`${appBase}/encyclopedia/#/search`}
        tooltip="Visit the NOMAD Materials Encyclopedia"
        icon={<EncyclopediaIcon/>}
      />}
    </MenuBarMenu>
    <MenuBarMenu name="analyze" route="/metainfo" icon={<AnalyticsIcon/>}>
      {!oasis && aitoolkitEnabled && <MenuBarItem
        label="AI Toolkit" name="aitoolkit" route="/aitoolkit"
        tooltip="NOMAD's Artificial Intelligence Toolkit tutorial Jupyter notebooks"
        icon={<MetainfoIcon />}
      />}
      <MenuBarItem
        name="metainfo" route="/metainfo" tooltip="Browse the NOMAD Archive schema"
      />
    </MenuBarMenu>
    <MenuBarMenu name="about" route="/" icon={<AboutIcon/>}>
      <MenuBarItem
        label="Information" name="about" route="/"
        tooltip="About the NOMAD Repository and Archive"
      />
      <MenuBarItem
        name="forum"
        href="https://matsci.org/c/nomad/"
        tooltip="The NOMAD user/developer forum on matsci.org"
        icon={<ForumIcon/>}
      />
      <MenuBarItem
        label="FAQ" name="faq"
        href="https://nomad-lab.eu/repository-archive-faqs"
        tooltip="Frequently Asked Questions (FAQ)"
        icon={<FAQIcon/>}
      />
      <MenuBarItem
        name="Docs"
        href={`${appBase}/docs/index.html`}
        tooltip="The full user and developer documentation"
        icon={<DocIcon/>}
      />
      <MenuBarItem
        name="Sources"
        href="https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR"
        tooltip="NOMAD's main Gitlab project"
        icon={<CodeIcon/>}
      />
      <Consent />
    </MenuBarMenu>
  </MenuBar>
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
    helpButton: {
      marginLeft: theme.spacing(1)
    },
    toolbar: {
      paddingRight: theme.spacing(3)
    },
    logo: {
      height: theme.spacing(7),
      marginRight: theme.spacing(2)
    },
    content: {
      marginTop: theme.spacing(14),
      flexGrow: 1,
      backgroundColor: theme.palette.background.default,
      width: '100%',
      overflow: 'auto'
    },
    barActions: {
      display: 'flex',
      alignItems: 'center'
    },
    barButton: {
      borderColor: theme.palette.getContrastText(theme.palette.primary.main),
      marginRight: 0
    },
    mainMenu: {
      marginLeft: theme.spacing(1)
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
    '/dataset': 'Dataset',
    '/aitoolkit': 'Artificial Intelligence Toolkit'
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
                  <MuiLink href="https://nomad-lab.eu">
                    <img alt="The NOMAD logo" className={classes.logo} src={`${guiBase}/nomad.png`}></img>
                  </MuiLink>
                  <Typography variant="h6" color="inherit" noWrap>
                    {selected(toolbarTitles)}
                  </Typography>
                  {help ? <HelpDialog color="inherit" maxWidth="md" classes={{root: classes.helpButton}} {...help}/> : ''}
                </div>
                <div className={classes.barActions}>
                  <LoginLogout color="primary" classes={{button: classes.barButton}} />
                </div>
              </Toolbar>
              <div className={classes.mainMenu} >
                <MainMenu />
              </div>
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
  'faq': {
    exact: true,
    path: '/faq',
    component: FAQ
  },
  'search': {
    exact: true,
    path: '/search',
    navPath: 'explore/search',
    component: SearchPage
  },
  'userdata': {
    exact: true,
    path: '/userdata',
    navPath: 'publish/userdata',
    component: UserdataPage
  },
  'entry': {
    path: '/entry/id',
    defaultNavPath: 'explore/search',
    component: EntryPage
  },
  'entry_query': {
    exact: true,
    path: '/entry/query',
    defaultNavPath: 'explore/search',
    component: EntryQuery
  },
  'entry_pid': {
    path: '/entry/pid',
    defaultNavPath: 'explore/search',
    component: ResolvePID
  },
  'dataset': {
    path: '/dataset/id',
    defaultNavPath: 'explore/search',
    component: DatasetPage
  },
  'dataset_doi': {
    path: '/dataset/doi',
    defaultNavPath: 'explore/search',
    component: ResolveDOI
  },
  'uploads': {
    exact: true,
    path: '/uploads',
    navPath: 'publish/uploads',
    component: UploadPage
  },
  'metainfo': {
    path: '/metainfo',
    navPath: 'analyze/metainfo',
    component: MetainfoPage
  },
  'aitoolkit': {
    path: '/aitoolkit',
    navPath: 'analyze/aitoolkit',
    component: AIToolkitPage
  },
  'about': {
    exact: true,
    path: '/',
    navPath: 'about/about',
    component: About
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
                  children={props => props.match && <route.component {...props} />}
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
