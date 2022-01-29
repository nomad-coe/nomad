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

import React, { useContext, useRef, useLayoutEffect, useMemo, useState, useCallback, createRef } from 'react'
import PropTypes from 'prop-types'
import { makeStyles, Card, CardContent, Box, Typography, Grid, Chip, Tooltip } from '@material-ui/core'
import grey from '@material-ui/core/colors/grey'
import ArrowRightIcon from '@material-ui/icons/ArrowRight'
import classNames from 'classnames'
import { useLocation, useRouteMatch, Link } from 'react-router-dom'
import { ErrorHandler } from '../ErrorHandler'

function escapeBadPathChars(s) {
  return s.replaceAll('$', '$0').replaceAll('?', '$1').replaceAll('#', '$2').replaceAll('%', '$3').replaceAll('\\', '$4')
}

function unescapeBadPathChars(s) {
  return s.replaceAll('$4', '\\').replaceAll('$3', '%').replaceAll('$2', '#').replaceAll('$1', '?').replaceAll('$0', '$')
}

export function formatSubSectionName(name) {
  // return name.startsWith('section_') ? name.slice(8) : name
  return name
}

export class Adaptor {
  constructor(e) {
    this.e = e

    if (new.target === Adaptor) {
      throw new TypeError('Cannot construct Abstract instances directly')
    }
  }

  isLoaded() {
    // Return false to signal to the Browser component that this adaptor needs to load some
    // data before we can call itemAdaptor. Adaptors that need to load data should call
    // lane.update() when they are done, to trigger rendering.
    return true
  }

  itemAdaptor(key) {
    // Gets the adaptor for the next lane. Will only be called if isLoaded returns true.
    return null
  }

  render() {
    return ''
  }
}

const useBrowserStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    flexFlow: 'column',
    margin: `-${theme.spacing(2)}px`,
    marginBottom: `-${theme.spacing(3)}px`
  },
  lanesContainer: {
    flex: '1 1 auto',
    height: '100%',
    overflowX: 'auto',
    overflowY: 'hidden',
    scrollBehavior: 'smooth'
  },
  lanes: {
    display: 'flex',
    overflow: 'hidden',
    height: '100%',
    overflowY: 'hidden',
    width: 'fit-content'
  }
}))
export const Browser = React.memo(function Browser({adaptor, form}) {
  const classes = useBrowserStyles()
  const rootRef = useRef()
  const outerRef = useRef()
  const innerRef = useRef()
  const { pathname } = useLocation()
  const { url } = useRouteMatch()

  useLayoutEffect(() => {
    function update() {
      const height = window.innerHeight - outerRef.current.getBoundingClientRect().top - 24
      rootRef.current.style.height = `${height}px`
      const scrollAmmount = innerRef.current.clientWidth - outerRef.current.clientWidth
      outerRef.current.scrollLeft = Math.max(scrollAmmount, 0)
    }
    if (url !== undefined) {
      update()
      window.addEventListener('resize', update)
      return () => window.removeEventListener('resize', update)
    }
  })

  const [, setRender] = useState(0)
  const update = useCallback(() => {
    setRender(current => current + 1)
  }, [setRender])

  const lanes = useRef([])

  if (url === undefined) {
    // Can happen when navigating to another tab, possibly with the browser's back/forward buttons
    // We want to keep the cached lanes, in case the user goes back to this tab, so return immediately.
    return
  }

  // Update the lanes
  const oldLanesByPath = {}
  lanes.current.forEach(lane => { oldLanesByPath[lane.path] = lane })

  const rootPath = url.endsWith('/') ? url.substring(0, url.length - 1) : url
  const root = oldLanesByPath[rootPath] || {
    key: 'root',
    path: rootPath,
    adaptor: adaptor,
    next: null,
    update: update
  }

  lanes.current = [root]
  root.next = null
  const segments = pathname.substring(url.length).split('/').filter(segment => segment)
  segments.forEach(segment => {
    const prev = lanes.current[lanes.current.length - 1]
    if (prev.adaptor) {
      const path = prev.path + '/' + encodeURI(segment)
      segment = unescapeBadPathChars(segment)
      const curr = oldLanesByPath[path] || {
        key: segment,
        path: path,
        data: root.adaptor.e,
        update: update
      }
      curr.next = undefined
      prev.next = curr
      if (!curr.adaptor && prev.adaptor.isLoaded()) {
        curr.adaptor = prev.adaptor.itemAdaptor(segment)
      }
      lanes.current.push(curr)
    }
  })

  return <React.Fragment>
    {form}
    <Card>
      <CardContent>
        <div className={classes.root} ref={rootRef} >
          <div className={classes.lanesContainer} ref={outerRef} >
            <div className={classes.lanes} ref={innerRef} >
              {lanes.current.map((lane, index) => (
                <Lane key={index} lane={lane} />
              ))}
            </div>
          </div>
        </div>
      </CardContent>
    </Card>
  </React.Fragment>
})
Browser.propTypes = ({
  adaptor: PropTypes.object.isRequired,
  form: PropTypes.node
})
export default Browser

export const laneContext = React.createContext()
const useLaneStyles = makeStyles(theme => ({
  root: {
    width: 'min-content',
    borderRight: `solid 1px ${grey[500]}`,
    display: 'inline-block'
  },
  container: {
    minWidth: 300,
    display: 'inline-block',
    height: '100%',
    overflowY: 'scroll'
  },
  error: {
    margin: theme.spacing(1)
  }
}))
function Lane({lane}) {
  const classes = useLaneStyles()
  const containerRef = createRef()
  const { key, adaptor, next } = lane
  lane.containerRef = containerRef
  const content = useMemo(() => {
    if (!adaptor) {
      return ''
    }
    return <div className={classes.root}>
      <div className={classes.container} ref={containerRef}>
        <laneContext.Provider value={lane}>
          <ErrorHandler message='This section could not be rendered, due to an unexpected error.' className={classes.error}>
            {adaptor.render()}
          </ErrorHandler>
        </laneContext.Provider>
      </div>
    </div>
    // We deliberetly break the React rules here. The goal is to only update if the
    // lanes contents change and not the lane object.
    // eslint-disable-next-line
  }, [key, adaptor, next?.key, classes])
  return content
}
Lane.propTypes = ({
  lane: PropTypes.object.isRequired
})

const useItemStyles = makeStyles(theme => ({
  root: {
    color: theme.palette.text.primary,
    textDecoration: 'none',
    margin: `0 -${theme.spacing(1)}px`,
    padding: `0 0 0 ${theme.spacing(1)}px`,
    whiteSpace: 'nowrap',
    display: 'flex',
    '& $icon': {
      color: theme.palette.grey[700]
    }
  },
  rootSelected: {
    backgroundColor: theme.palette.primary.main,
    color: theme.palette.primary.contrastText,
    whiteSpace: 'nowrap',
    '& $icon': {
      color: theme.palette.primary.contrastText
    }
  },
  rootUnSelected: {
    '&:hover': {
      backgroundColor: grey[300]
    }
  },
  disabled: {
    color: 'grey'
  },
  childContainer: {
    flexGrow: 1,
    overflow: 'hidden',
    textOverflow: 'ellipsis',
    whiteSpace: 'nowrap'
  },
  icon: {
    fontSize: 16,
    marginTop: 4
  }
}))

export function Item({children, itemKey, disabled, icon, actions, chip}) {
  const classes = useItemStyles()
  const lane = useContext(laneContext)
  const selected = lane.next && lane.next.key
  const isSelected = selected === itemKey
  if (disabled) {
    return <div className={classNames(classes.childContainer, classes.disabled)}>{children}</div>
  }
  return <Link
    className={classNames(
      classes.root,
      isSelected ? classes.rootSelected : classes.rootUnSelected
    )}
    to={`${lane.path}/${encodeURI(escapeBadPathChars(itemKey))}`}
  >
    <Grid container spacing={2} alignItems="center">
      {icon && <Grid item>
        {React.createElement(icon, {fontSize: 'small', className: classes.icon})}
      </Grid>}
      <Grid item className={classes.childContainer}>
        <Typography>{children}</Typography>
      </Grid>
      {chip && (
        <Grid item>
          <Chip size="small" color={isSelected ? 'primary' : 'default'} label={chip} />
        </Grid>
      )}
      {actions && <Grid item>
        {actions}
      </Grid>}
    </Grid>
    <ArrowRightIcon/>
  </Link>
}
Item.propTypes = ({
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]).isRequired,
  itemKey: PropTypes.string.isRequired,
  disabled: PropTypes.bool,
  icon: PropTypes.elementType,
  chip: PropTypes.string,
  actions: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
})

export function Content(props) {
  return <Box maxWidth={1024} padding={1} {...props} />
}

export function Compartment({title, children, color}) {
  if (!React.Children.count(children)) {
    return ''
  }
  return <React.Fragment>
    <Box paddingTop={1} whiteSpace="nowrap">
      {title && <Typography color={color} variant="overline">{title}</Typography>}
    </Box>
    {children}
  </React.Fragment>
}
Compartment.propTypes = ({
  title: PropTypes.string,
  color: PropTypes.string,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]).isRequired
})

export function Title({title, label, tooltip, actions, ...moreProps}) {
  return <Compartment>
    <Grid container justifyContent="space-between" wrap="nowrap" spacing={1}>
      <Grid item>
        {tooltip ? (
          <Tooltip title={tooltip}>
            <div style={{overflow: 'hidden', textOverflow: 'ellipsis'}}>
              <Typography variant="h6" {...moreProps}>{title}</Typography>
            </div>
          </Tooltip>
        ) : (
          <Typography variant="h6" {...moreProps}>{title}</Typography>
        )}
        {label && (
          <Typography variant="caption" color={moreProps.color}>
            {label}
          </Typography>
        )}
      </Grid>
      {actions && (
        <Grid item>
          {actions}
        </Grid>
      )}
    </Grid>
  </Compartment>
}
Title.propTypes = ({
  title: PropTypes.string.isRequired,
  label: PropTypes.string,
  tooltip: PropTypes.string,
  actions: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
})
