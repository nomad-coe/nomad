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

import React, { useCallback, useMemo } from 'react'
import { matchPath, useLocation, Link as RouterLink } from 'react-router-dom'
import { Typography, Breadcrumbs as MUIBreadcrumbs, Link, Box, makeStyles, Tooltip } from '@material-ui/core'
import { HelpButton } from '../Help'
import { allRoutes } from './Routes'
import {useDataStore} from "../DataStore"

const useStyles = makeStyles(theme => ({
  root: {
    marginLeft: 5,
    marginTop: theme.spacing(-0.25)
  },
  help: {
    marginLeft: theme.spacing(0.5)
  },
  ellipsis: {
    direction: 'rtl',
    textAlign: 'left',
    whiteSpace: 'nowrap',
    overflow: 'hidden',
    textOverflow: 'ellipsis',
    maxWidth: 350
  }
}))

const Breadcrumbs = React.memo(function Breadcrumbs() {
  const styles = useStyles()
  const dataStore = useDataStore()
  const {pathname} = useLocation()
  const routes = useMemo(() => allRoutes.slice().reverse(), [])

  const getRef = useCallback((breadcrumb) => {
    if (breadcrumb === 'Entry') {
      return dataStore.breadcrumb.entryRef
    } else if (breadcrumb === 'Upload') {
      return dataStore.breadcrumb.uploadRef
    } else {
      return undefined
    }
  }, [dataStore.breadcrumb.entryRef, dataStore.breadcrumb.uploadRef])

  const getFinalBreadcrumb = useCallback((breadcrumb) => {
    if (breadcrumb === 'Entry') {
      return dataStore.breadcrumb.getEntry()
    } else if (breadcrumb === 'Upload') {
      return dataStore.breadcrumb.getUpload()
    } else {
      return breadcrumb
    }
  }, [dataStore.breadcrumb])

  return <MUIBreadcrumbs className={styles.root}>
    {routes
      .filter(route => route.breadcrumb && route.path)
      .map(route => ({route: route, match: matchPath(pathname, {path: route.path})}))
      .filter(({match}) => match)
      .map(({route, match}, i) => {
        const title = <Typography
          className={styles.ellipsis}
          variant={i ? undefined : 'h6'}
          color="textPrimary"
          ref={getRef(route.breadcrumb)}
        >
          {getFinalBreadcrumb(route.breadcrumb)}
        </Typography>
        if (match.url === pathname) {
          return <Box key={i} display="flex" flexDirection="row" alignItems="center">
            {title}
            {route.help && (
              <Tooltip title={route?.help?.title || ""}>
                <HelpButton className={styles.help} size="small" IconProps={{fontSize: 'small'}} heading={route?.help?.title} content={route?.help?.content} />
              </Tooltip>
            )}
          </Box>
        } else {
          return <Link key={i} component={RouterLink} to={match.url}>{title}</Link>
        }
      })}
  </MUIBreadcrumbs>
})

export default Breadcrumbs
