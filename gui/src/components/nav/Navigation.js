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

import React, { useEffect, useRef, useState } from 'react'
import { makeStyles } from '@material-ui/core/styles'
import { Snackbar, SnackbarContent, IconButton, Link as MuiLink, Button } from '@material-ui/core'
import UnderstoodIcon from '@material-ui/icons/Check'
import ReloadIcon from '@material-ui/icons/Replay'
import { amber } from '@material-ui/core/colors'
import AppBar from './AppBar'
import { guiBase, version } from '../../config'
import packageJson from '../../../package.json'
import Routes from './Routes'
import { withApi } from '../api'
import { reload } from '../../serviceWorker'

export const ScrollContext = React.createContext({scrollParentRef: null})

function ReloadSnack() {
  return <Snackbar
    anchorOrigin={{
      vertical: 'bottom',
      horizontal: 'left'
    }}
    open
  >
    <SnackbarContent
      message={<span>There is a new NOMAD version. Please reload the app.</span>}
      action={[
        <Button
          key={0} color="inherit" startIcon={<ReloadIcon/>}
          onClick={() => reload()}
        >
          reload
        </Button>
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
  const [showReloadSnack, setShowReloadSnack] = useState(false)
  const scrollParentRef = useRef(null)

  useEffect(() => {
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
          console.log('GUI version mismatch')
          setShowReloadSnack(true)
        }
      })
      .catch(() => {
        console.log('Could not validate version, continue...')
      })
  }, [setShowReloadSnack])

  return (
    <div className={classes.root}>
      <div className={classes.appFrame}>
        <BetaSnack />
        { showReloadSnack ? <ReloadSnack/> : ''}
        <AppBar />

        <main className={classes.content} ref={scrollParentRef}>
          <ScrollContext.Provider value={{scrollParentRef: scrollParentRef}}>
            <Routes/>
          </ScrollContext.Provider>
        </main>
      </div>
    </div>
  )
}

export default withApi(false)(Navigation)
