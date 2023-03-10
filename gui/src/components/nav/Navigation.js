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
import AppBar, { appBarHeight } from './AppBar'
import { guiBase, version } from '../../config'
import { Routes } from './Routes'
import { usePWA } from '../PWA'
import { ErrorBoundary } from '../errors'
import { useCookies } from 'react-cookie'

export const ScrollContext = React.createContext({scrollParentRef: null})

function ReloadSnack() {
  const {showReload, reloadPage} = usePWA()

  return <Snackbar
    anchorOrigin={{
      vertical: 'bottom',
      horizontal: 'left'
    }}
    open={showReload}
  >
    <SnackbarContent
      message={<span>There is a new NOMAD version. Please reload the app.</span>}
      action={[
        <Button
          key={0} color="inherit" startIcon={<ReloadIcon/>}
          onClick={() => reloadPage()}
        >
          reload
        </Button>
      ]}
    />
  </Snackbar>
}

const useTermsSnackStyles = makeStyles(theme => ({
  termsLink: {
    color: theme.palette.primary.light
  }
}))

function TermsSnack() {
  const [cookies, setCookie] = useCookies()
  const [accepted, setAccepted] = useState(cookies['terms-accepted'])
  const classes = useTermsSnackStyles()

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
        our <Link className={classes.termsLink} href="https://nomad-lab.eu/nomad-lab/terms.html" title="terms of service">terms of service</Link>.
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
    return null
  }

  if (!version.isBeta && !version.isTest) {
    return null
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
    marginTop: theme.spacing(appBarHeight),
    flexGrow: 1,
    backgroundColor: theme.palette.background.default,
    width: '100%',
    overflow: 'auto'
  }
}))

export default function Navigation() {
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
