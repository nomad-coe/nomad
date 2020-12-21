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

import React, { useCallback, useContext, useEffect, useState } from 'react'
import { useLocation } from 'react-router-dom'
import { Typography,
  AppBar as MuiAppBar, Toolbar, Link, LinearProgress, makeStyles } from '@material-ui/core'
import { allRoutes as routes } from './Routes'
import LoginLogout from '../LoginLogout'
import HelpDialog from '../Help'
import MainMenu from './MainMenu'
import { apiContext } from '../api'
import { guiBase } from '../../config'

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

const useAppBarStyles = makeStyles(theme => ({
  root: {
    zIndex: theme.zIndex.drawer + 1,
    backgroundColor: 'white'
  },
  title: {
    marginLeft: theme.spacing(1),
    flexGrow: 1,
    display: 'flex',
    alignItems: 'center',
    alignContent: 'flex-start',
    color: theme.palette.primary.main
  },
  toolbar: {
    paddingRight: theme.spacing(3)
  },
  helpButton: {
    marginLeft: theme.spacing(1)
  },
  logo: {
    height: theme.spacing(7),
    marginRight: theme.spacing(2)
  },
  actions: {
    display: 'flex',
    alignItems: 'center'
  },
  button: {
    borderColor: theme.palette.getContrastText(theme.palette.primary.main),
    marginRight: 0
  },
  mainMenu: {
    marginLeft: theme.spacing(1)
  }
}))

export default function AppBar() {
  const classes = useAppBarStyles()
  const {pathname} = useLocation()
  const selectedRouteKey = Object.keys(routes).find(key => {
    const route = routes[key]
    return pathname.startsWith(route.path)
  })
  const selectedRoute = routes[selectedRouteKey]
  const help = selectedRoute.appBarHelp
  const title = selectedRoute.appBarTitle

  return <MuiAppBar position="fixed" className={classes.root}>
    <Toolbar classes={{root: classes.toolbar}} disableGutters>
      <div className={classes.title}>
        <Link href="https://nomad-lab.eu">
          <img alt="The NOMAD logo" className={classes.logo} src={`${guiBase}/nomad.png`}></img>
        </Link>
        <Typography variant="h6" color="inherit" noWrap>
          {title}
        </Typography>
        {help ? <HelpDialog color="inherit" maxWidth="md" classes={{root: classes.helpButton}} {...help}/> : ''}
      </div>
      <div className={classes.actions}>
        <LoginLogout color="primary" classes={{button: classes.button}} />
      </div>
    </Toolbar>
    <div className={classes.mainMenu} >
      <MainMenu />
    </div>
    <LoadingIndicator />
  </MuiAppBar>
}
