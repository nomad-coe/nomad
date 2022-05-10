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

import React, { useEffect, useMemo, useRef, useState } from 'react'
import { makeStyles } from '@material-ui/core/styles'
import { useLocation } from 'react-router-dom'
import { Snackbar, SnackbarContent, IconButton, Link as MuiLink, Button, Link } from '@material-ui/core'
import UnderstoodIcon from '@material-ui/icons/Check'
import ReloadIcon from '@material-ui/icons/Replay'
import { amber } from '@material-ui/core/colors'
import AppBar from './AppBar'
import { version, guiBase } from '../../config'
import Routes from './Routes'
import { withApi } from '../api'
import { serviceWorkerUpdateHandlerRef } from '../../index'
import { ErrorBoundary } from '../errors'
import { useCookies } from 'react-cookie'

export const ScrollContext = React.createContext({scrollParentRef: null})

function ReloadSnack() {
  const waitingServiceWorker = useRef(null)
  const [reload, setReload] = useState(false)
  useEffect(() => {
    serviceWorkerUpdateHandlerRef.current = registration => {
      waitingServiceWorker.current = registration.waiting
      setReload(true)
    }
  }, [setReload, waitingServiceWorker])

  return <Snackbar
    anchorOrigin={{
      vertical: 'bottom',
      horizontal: 'left'
    }}
    open={reload}
  >
    <SnackbarContent
      message={<span>There is a new NOMAD version. Please reload the app.</span>}
      action={[
        <Button
          key={0} color="inherit" startIcon={<ReloadIcon/>}
          onClick={() => {
            if (waitingServiceWorker.current) {
              waitingServiceWorker.current.onstatechange = () => {
                if (waitingServiceWorker.current.state === 'activated') {
                  window.location.reload()
                }
              }
              waitingServiceWorker.current.postMessage({type: 'SKIP_WAITING'})
            }
          }}
        >
          reload
        </Button>
      ]}
    />
  </Snackbar>
}

function TermsSnack() {
  const [cookies, setCookie] = useCookies()
  const [accepted, setAccepted] = useState(cookies['terms-accepted'])

  const cookieOptions = useMemo(() => ({
    expires: new Date(2147483647 * 1000),
    path: '/' + guiBase.split('/').slice(1).join('/')
  }), [])

  return <Snackbar
    anchorOrigin={{
      vertical: 'bottom',
      horizontal: 'left'
    }}
    open={!accepted}
  >
    <SnackbarContent
      message={<span>
        NOMAD only uses cookies that are strictly necessary for this site&apos;s functionality.
        No tracking or marketing cookies are used. By using this site you agree to
        our <Link href="https://nomad-lab.eu/terms" title="terms of service">terms of service</Link>.
      </span>}
      action={[
        <IconButton
          size="small" key={0} color="inherit"
          onClick={() => {
            setCookie('terms-accepted', true, cookieOptions)
            setAccepted(true)
          }}
        >
          <UnderstoodIcon />
        </IconButton>
      ]}
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
        <IconButton size="small" key={0} color="inherit" onClick={() => setUnderstood(true)}>
          <UnderstoodIcon />
        </IconButton>
      ]}
    />
  </Snackbar>
}

function LegacySnack() {
  const classes = useBetaSnackStyles()
  const [understood, setUnderstood] = useState(false)

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
        You are using an old version of NOMAD. This version can still be used to view data.
        To change yor data or upload new data, please use
        <MuiLink style={{color: 'white'}} href="https://nomad-lab.eu/prod/v1/gui">the official NOMAD version</MuiLink>.
      </span>}
      action={[
        <IconButton size="small" key={0} color="inherit" onClick={() => setUnderstood(true)}>
          <UnderstoodIcon />
        </IconButton>
      ]}
    />
  </Snackbar>
}

const useStyles = makeStyles(theme => ({
  root: {
    minWidth: 1024
  },
  appFrame: {
    zIndex: 1,
    overflow: 'hidden',
    position: 'relative',
    display: 'flex',
    width: '100%',
    height: '100vh'
  },
  content: {
    marginTop: theme.spacing(14),
    flexGrow: 1,
    backgroundColor: theme.palette.background.default,
    width: '100%',
    overflow: 'auto'
  }
}))

function Navigation() {
  const classes = useStyles()
  const { pathname } = useLocation()
  const scrollParentRef = useRef(null)

  // Scroll to top upon changing page
  useEffect(() => {
    if (scrollParentRef.current) {
      scrollParentRef.current.scrollTo(0, 0)
    }
  }, [pathname])

  return (
    <div className={classes.root}>
      <div className={classes.appFrame}>
        <ReloadSnack/>
        <ErrorBoundary>
          <LegacySnack />
          <BetaSnack />
          <TermsSnack />
          <AppBar />
          <main className={classes.content} ref={scrollParentRef}>
            <ScrollContext.Provider value={{scrollParentRef: scrollParentRef}}>
              <Routes/>
            </ScrollContext.Provider>
          </main>
        </ErrorBoundary>
      </div>
    </div>
  )
}

export default withApi(false)(Navigation)
