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

import React, {useCallback, useContext, useEffect, useLayoutEffect, useMemo, useRef, useState} from 'react'
import PropTypes from 'prop-types'
import { makeStyles, Card, CardContent, Box, Typography, Grid, Chip, Tooltip, IconButton } from '@material-ui/core'
import grey from '@material-ui/core/colors/grey'
import classNames from 'classnames'
import { useLocation, useRouteMatch, Link } from 'react-router-dom'
import { ErrorHandler } from '../ErrorHandler'
import { useDataStore } from '../DataStore'
import { useApi } from '../api'
import NavigateIcon from '@material-ui/icons/ArrowRight'
import H5Web from '../visualization/H5Web'
import {CopyToClipboard} from "react-copy-to-clipboard"
import ClipboardIcon from "@material-ui/icons/Assignment"

function escapeBadPathChars(s) {
  return s.replace(/!/g, '!0').replace(/\?/g, '!1').replace(/#/g, '!2').replace(/%/g, '!3').replace(/\\/g, '!4')
}

function unescapeBadPathChars(s) {
  return s.replace(/!4/g, '\\').replace(/!3/g, '%').replace(/!2/g, '#').replace(/!1/g, '?').replace(/!0/g, '!')
}

/**
 * Browsers are made out of lanes. Each lane uses an adaptor that determines how to render
 * the lane contents and what adaptor is used for the next lane (depending on what is
 * selected in this lane).
 */
export class Adaptor {
  constructor() {
    if (new.target === Adaptor) {
      throw new TypeError('Cannot construct Abstract instances directly')
    }
  }

  /**
   * A potentially asynchronous method called to initialize an adaptor, before any
   * calls are made to the itemAdaptor or render methods.
   */
  initialize(api, dataStore) {
  }

  /**
   * A potentially asynchronous method called when the adaptor is cleaned up due to user navigation
   */
  cleanup() {
  }

  /**
   * Optional additional methods:
   *   onScrollToEnd: called when the user scrolls to the bottom of the lane. Useful to load more info.
   *   onRendered: called in a useEffect after the adaptor has been rendered.
   */

  /**
   * A potentially asynchronous method that is used to determine the adaptor for the
   * next lane depending on the given key/url segment.
   * @returns An adaptor that is used to render the next lane.
   */
  itemAdaptor(key) {
    return null
  }

  /**
   * Optionally returns a set of strings defining this lanes dependencies. These dependencies
   * are used to determine if an adaptor is to be invalidated when calling invalidateLanesWithDependency
   */
  depends() {
    return null
  }

  /**
   * Renders the contents of the current lane (the lane that this adaptor represents).
   */
  render() {
    return null
  }
}

export const browserContext = React.createContext()

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
  },
  sideLane: {
    display: 'flex',
    height: '100%'
  }
}))
export const Browser = React.memo(function Browser({adaptor, form}) {
  const classes = useBrowserStyles()
  const dataStore = useDataStore()
  const rootRef = useRef()
  const outerRef = useRef()
  const innerRef = useRef()
  const { pathname } = useLocation()
  const { url } = useRouteMatch()

  const { api } = useApi()

  const [hdf5Path, setHdf5Path] = useState(null)
  const [hdf5Filename, setHdf5Filename] = useState(null)

  const checkHdf5Path = (m_nx_data_path) => {
    if (typeof m_nx_data_path !== 'undefined') {
      setHdf5Path(m_nx_data_path)
    } else {
      setHdf5Path(null)
    }
  }

  const checkHdf5File = (m_nx_data_file) => {
    if (typeof m_nx_data_file !== 'undefined') {
      setHdf5Filename(m_nx_data_file)
    }
  }

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

  const [, setInternalRender] = useState(0)
  const internalUpdate = useCallback(() => {
    setInternalRender(current => current + 1)
  }, [setInternalRender])
  const lanes = useRef([])
  const computingForUrl = useRef(null)
  const computeLanes = useCallback(async () => {
    if (!url || computingForUrl.current === url) {
      // We have no url (happens if we navigate to a different tab), or we're already in the
      // midst of computing lanes for this url
      return
    }
    computingForUrl.current = url
    const rootPath = url.endsWith('/') ? url.substring(0, url.length - 1) : url
    const path = pathname?.replace(/\/(-?\d*)(\/|$)/g, ":$1/").replace(/\/$/, "")
    const segments = ['root'].concat(path.substring(url.length).split('/').filter(segment => segment))
    const oldLanes = lanes.current
    const newLanes = []

    for (let index = 0; index < segments.length; index++) {
      const segment = unescapeBadPathChars(segments[index])
      const prev = index === 0 ? null : newLanes[index - 1]
      let lane
      let newLaneCreated = false

      if (oldLanes[index]?.key === segment && !newLaneCreated) {
        // reuse the existing lane (incl. its adaptor and data)
        lane = oldLanes[index]
        lane.next = null
      } else {
        // create new lane
        lane = {
          index: index,
          key: segment,
          path: prev ? prev.path + '/' + encodeURI(escapeBadPathChars(segment)) : rootPath,
          next: null,
          prev: prev,
          initialized: false
        }
        if (!newLaneCreated) {
          // First newly created lane. Cleanup all adaptors that are being discarded
          oldLanes.slice(index).forEach(oldLane => oldLane.adaptor?.cleanup())
        }
        newLaneCreated = true

        if (typeof adaptor?.obj?.metadata?.mainfile !== 'undefined') {
          setHdf5Filename(adaptor?.obj?.metadata?.mainfile)
        }

        // Set or create the lane adaptor
        if (!prev) {
          // The root lane - use provided adaptor
          lane.adaptor = adaptor
        } else {
          // Non-root lane - create adaptor
          try {
            lane.adaptor = await prev.adaptor.itemAdaptor(segment)
            checkHdf5Path(lane.adaptor?.obj?.m_attributes?.m_nx_data_path)
            checkHdf5File(lane.adaptor?.obj?.m_attributes?.m_nx_data_file)
          } catch (error) {
            console.log(error)
            lane.error = `The item "${segment}" could not be found.`
          }
        }

        // initialize the adaptor
        if (lane.adaptor) {
          try {
            await lane.adaptor.initialize(api, dataStore)
            lane.initialized = true
          } catch (error) {
            console.log(error)
            lane.error = 'Could not initialize view' + (lane.index > 0 ? `: item "${segment}" not found.` : '.')
          }
        }
      }

      if (prev) {
        prev.next = lane
        checkHdf5Path(lane.adaptor?.obj?.m_attributes?.m_nx_data_path)
        checkHdf5File(lane.adaptor?.obj?.m_attributes?.m_nx_data_file)
      }
      newLanes.push(lane)
      if (lane.error) {
        break // Ignore subsequent segments/lanes
      }
    }
    lanes.current = newLanes
    computingForUrl.current = null
    internalUpdate()
  }, [adaptor, api, dataStore, internalUpdate, url, pathname])

  // Method used to invalidate and refresh lanes from the provided lane index and forward.
  const invalidateLanesFromIndex = useCallback((index) => {
    index = index || 0
    const droppedLanes = lanes.current.splice(index)
    droppedLanes.forEach(droppedLane => droppedLane.adaptor?.cleanup())
    computeLanes()
  }, [computeLanes])

  // Method used to invalidate all lanes that have a certain dependency. The first lane which
  // has this dependency, and all subsequent lanes, are invalidated.
  const invalidateLanesWithDependency = useCallback((dependency) => {
    for (const lane of lanes.current) {
      const dependencies = lane.adaptor && lane.adaptor.depends()
      if (dependencies && dependencies.has(dependency)) {
        invalidateLanesFromIndex(lane.index)
        return
      }
    }
  }, [invalidateLanesFromIndex])

  useEffect(() => {
    return () => {
      // Cleanup method, will run when the root adaptor changes or the Browser component is dismounted
      // We should then cleanup all existing adaptors and not reuse any lanes.
      lanes.current.forEach(lane => lane.adaptor?.cleanup())
      lanes.current = []
    }
  }, [adaptor])

  useEffect(() => {
    computeLanes()
  }, [computeLanes])

  const contextValue = useMemo(() => {
    return {
      lanes,
      invalidateLanesFromIndex,
      invalidateLanesWithDependency
    }
  }, [lanes, invalidateLanesFromIndex, invalidateLanesWithDependency])

  if (url === undefined) {
    // Can happen when navigating to another tab, possibly with the browser's back/forward buttons
    // We want to keep the cached lanes, in case the user goes back to this tab, so return immediately.
    return null
  }

  return <browserContext.Provider value={contextValue}>
    {form}
      <Grid container direction="row" spacing={2}>
        <Grid item md={((hdf5Path && hdf5Filename) ? 8 : 12)}>
        <Card>
          <CardContent>
        <div className={classes.root} ref={rootRef}>
          <div className={classes.lanesContainer} ref={outerRef}>
            <div className={classes.lanes} ref={innerRef}>
              {lanes.current && lanes.current.map((lane, index) => (
                <Lane key={`lane${index}`} lane={lane} />
              ))}
            </div>
          </div>
        </div>
        </CardContent>
        </Card>
        </Grid>
        {hdf5Path && hdf5Filename && adaptor?.obj?.metadata?.upload_id && <Grid item md={4} style={{visibility: (hdf5Path == null ? "hidden" : "visible")}}>
          <Card className={classes.sideLane}>
            {/* We use the key prop in the H5Web component here to force re-render on path changes while browsing with the browser. The prop, initialPath, does not re-render the component. */}
            <H5Web key={hdf5Path} upload_id={adaptor?.obj?.metadata?.upload_id} filename={hdf5Filename} initialPath={hdf5Path} sidebarOpen={false}/>
          </Card>
        </Grid>}
      </Grid>
  </browserContext.Provider>
})
Browser.propTypes = ({
  adaptor: PropTypes.object.isRequired,
  form: PropTypes.node
})
export default Browser

export const laneContext = React.createContext()

export const useLane = () => {
  return useContext(laneContext)
}

export const laneErrorBoundryMessage = 'This section could not be rendered, due to an unexpected error'

const useLaneStyles = makeStyles(theme => ({
  root: {
    width: 'min-content',
    borderRight: `solid 1px ${grey[500]}`,
    display: 'inline-block'
  },
  container: {
    display: 'inline-block',
    height: '100%',
    overflowX: 'clip',
    overflowY: 'scroll'
  },
  error: {
    margin: theme.spacing(1),
    minWidth: 300
  }
}))
function Lane({lane}) {
  const classes = useLaneStyles()
  const containerRef = useRef(null)
  const { key, adaptor, next, initialized, error } = lane
  const oldScrollValue = useRef(0)
  lane.containerRef = containerRef
  const [internalUpdate, setInternalRender] = useState(0)
  const updateLane = useCallback(() => {
    setInternalRender(current => current + 1)
  }, [setInternalRender])

  const handleScroll = useCallback((e) => {
    if (adaptor.onScrollToEnd) {
      const threshold = (e.target.scrollHeight - e.target.clientHeight) * 95 / 100
      if (e.target.scrollTop >= threshold && oldScrollValue.current < threshold) {
        adaptor.onScrollToEnd(e, updateLane)
      }
      oldScrollValue.current = e.target.scrollTop
    }
  }, [adaptor, updateLane])

  const content = useMemo(() => {
    if (error) {
      return <div className={classes.error}>
        <Typography variant="h6" color="error">ERROR</Typography>
        <Typography color="error">{error}</Typography>
      </div>
    }
    if (!adaptor || !initialized) {
      return null
    }
    return <div className={classes.root} key={`lane:${lane.path}`} data-testid={`lane${lane.index}:${lane.key}`}>
      <div className={classes.container} ref={containerRef} onScroll={handleScroll}>
        <laneContext.Provider value={lane}>
          <ErrorHandler message={laneErrorBoundryMessage} className={classes.error}>
            {adaptor.render()}
          </ErrorHandler>
        </laneContext.Provider>
      </div>
    </div>
    // We deliberetly break the React rules here. The goal is to only update if the
    // lanes contents change and not the lane object.
    // eslint-disable-next-line
  }, [lane, key, adaptor, initialized, next?.key, classes, error, internalUpdate])

  useEffect(() => {
    if (adaptor && initialized && containerRef.current && adaptor.onRendered) {
      adaptor.onRendered({target: containerRef.current}, updateLane)
    }
    // eslint-disable-next-line
  }, [key, adaptor, internalUpdate, updateLane])

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
  rootUnSelectedHighlighted: {
    color: theme.palette.primary.main,
    '& $icon': {
      color: theme.palette.primary.main
    },
    '&:hover': {
      backgroundColor: grey[300]
    }
  },
  disabled: {
    color: 'grey'
  },
  rightPaddedItem: {
    padding: `0 ${theme.spacing(1)}px 0 0`
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

const ItemLink = React.forwardRef(function ItemLink({itemKey, ...props}, ref) {
  const lane = useContext(laneContext)
  if (lane && itemKey) {
    return <Link
      {...props}
      data-testid={`item:${itemKey}`}
      to={`${lane.path}/${encodeURI(escapeBadPathChars(itemKey))}`}
    />
  } else {
    // If this is used in a non Browser context
    return <div
      {...props}
      data-testid={`item:${itemKey}`}
    />
  }
})
ItemLink.propTypes = {
  itemKey: PropTypes.string.isRequired
}

export function ItemButton({itemKey, itemLink, icon, ...props}) {
  // ItemButtons are often used in clickable context and we should stop the propagation
  // of click events to prevent unwanted behavior.
  const handleClick = useCallback((event) => {
    event.stopPropagation()
  }, [])
  return (
      <IconButton {...props} component={itemLink || ItemLink} itemKey={itemKey} onClick={handleClick}>
        {icon || <NavigateIcon/>}
      </IconButton>
  )
}
ItemButton.propTypes = {
  itemKey: PropTypes.string.isRequired,
  itemLink: PropTypes.element,
  icon: PropTypes.node
}

export function Item({children, itemKey, length, disabled, highlighted, icon, actions, chip}) {
  const classes = useItemStyles()
  const lane = useLane()
  let selected = lane?.next && lane?.next.key
  let [label, index] = selected ? selected.split(':') : []
  if (index && length) {
    index = parseInt(index)
    if (index < 0) index = index + length
    selected = `${label}:${index}`
  }
  const isSelected = itemKey && (selected === itemKey || selected?.replace(':', '/') === itemKey)
  if (disabled) {
    return <Grid
      container spacing={0} alignItems="center" wrap="nowrap"
      style={{padding: 0, margin: 0}}
    >
      {icon && <Grid item className={classes.rightPaddedItem}>
        {React.createElement(icon, {fontSize: 'small', className: classes.icon})}
      </Grid>}
      <Grid item className={classNames(classes.childContainer, classes.disabled, classes.rightPaddedItem)}>
        <Typography noWrap>{children}</Typography>
      </Grid>
      {actions && <Grid item>
        {actions}
      </Grid>}
    </Grid>
  }

  return <ItemLink
    className={classNames(
      classes.root,
      isSelected ? classes.rootSelected : highlighted ? classes.rootUnSelectedHighlighted : classes.rootUnSelected
    )}
    itemKey={itemKey}
  >
    <Grid
      container spacing={0} alignItems="center" wrap="nowrap"
      style={{padding: 0, margin: 0}}
    >
      {icon && <Grid item className={classes.rightPaddedItem}>
        {React.createElement(icon, {fontSize: 'small', className: classes.icon})}
      </Grid>}
      <Grid item className={classNames(classes.childContainer, classes.rightPaddedItem)}>
        {children}
      </Grid>
      {chip && (
        <Grid item className={classes.rightPaddedItem}>
          <Chip size="small" color={isSelected ? 'primary' : 'default'} label={chip} />
        </Grid>
      )}
      {actions && <Grid item className={classes.rightPaddedItem}>
        {actions}
      </Grid>}
      {itemKey && <Grid item style={{padding: 0}}>
        <NavigateIcon padding="0"/>
      </Grid>}
    </Grid>
  </ItemLink>
}
Item.propTypes = ({
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]).isRequired,
  itemKey: PropTypes.string,
  length: PropTypes.number,
  disabled: PropTypes.bool,
  highlighted: PropTypes.bool,
  icon: PropTypes.elementType,
  chip: PropTypes.string,
  actions: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]),
  open: PropTypes.bool,
  onClick: PropTypes.func
})

export function Content(props) {
  return <Box minWidth={300} maxWidth={600} padding={1} {...props} />
}

export function Compartment({title, children, color, startCollapsed, onUnfold}) {
  const [collapsed, setCollapsed] = useState(startCollapsed)
  const handleClick = useCallback(() => {
    setCollapsed(value => !value)
    onUnfold && onUnfold()
  }, [setCollapsed, onUnfold])
  if (!React.Children.count(children)) {
    return null
  }

  return <React.Fragment>
    <Box paddingTop={1} whiteSpace="nowrap" onClick={(handleClick)} style={{cursor: 'pointer'}} data-testid={'compartment'}>
      {title && <Typography color={color} variant="overline">{title}</Typography>}
      {collapsed && <ItemChip label="closed" color="primary" data-testid={`collapsed:${title}`}/>}
    </Box>
    {(!collapsed) && children}
  </React.Fragment>
}
Compartment.propTypes = ({
  title: PropTypes.string,
  color: PropTypes.string,
  startCollapsed: PropTypes.bool,
  onUnfold: PropTypes.func,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ]).isRequired
})

const useTitleStyles = makeStyles(theme => ({
  titleGridItem: {
    flexGrow: 1,
    overflow: 'hidden',
    textOverflow: 'ellipsis'
  },
  sectionLabel: {
    whiteSpace: 'nowrap'
  },
  sectionValue: {
    whiteSpace: 'nowrap',
    overflow: 'hidden',
    display: 'flex',
    marginLeft: 3
  },
  sectionName: {
    overflow: 'hidden',
    whiteSpace: 'nowrap',
    textOverflow: 'ellipsis',
    flexShrink: 1
  }
}))

function LinkWithCopyToClipboard(props) {
  const classes = useTitleStyles()
  const {label, name, itemKey, ...moreProps} = props
  const lane = useLane()
  const selected = lane?.next && lane?.next.key
  const isSelected = itemKey && (selected === itemKey)

  return <Box display={'flex'} alignItems={'center'} marginBottom={1}>
    {label && (
      <Typography className={classes.sectionLabel} color={moreProps.color}>
        {label}
      </Typography>
    )}
    {name && (
      <Box className={classes.sectionValue}>
        <ItemLink itemKey={itemKey} className={classes.sectionValue} style={{ textDecoration: 'none' }}>
          {isSelected
            ? <Chip
              classes={{label: classes.sectionName}}
              style={{maxWidth: '100%'}}
              label={name}
              color="primary"
              size="small"
            />
            : <Typography className={classes.sectionName} color={'primary'}>
              {name}
            </Typography>}
        </ItemLink>
        <CopyToClipboard
          text={name}
          onCopy={() => null}
        >
          <Tooltip title={`Copy "${name}" to clipboard`}>
            <IconButton style={{ padding: 0, marginLeft: 3 }}>
              <ClipboardIcon fontSize="small"/>
            </IconButton>
          </Tooltip>
        </CopyToClipboard>
      </Box>
    )}
  </Box>
}
LinkWithCopyToClipboard.propTypes = ({
  label: PropTypes.string,
  name: PropTypes.string,
  itemKey: PropTypes.string
})

export function Title({title, label, definitionName, subSectionName, tooltip, actions, ...moreProps}) {
  const classes = useTitleStyles()
  return <Compartment>
    <Grid container justifyContent="space-between" wrap="nowrap" spacing={1}>
      <Grid item className={classes.titleGridItem}>
        {tooltip ? (
          <Tooltip title={tooltip}>
            <Typography variant="h6" noWrap {...moreProps}>{title || <i>unnamed</i>}</Typography>
          </Tooltip>
        ) : (
          <Typography variant="h6" noWrap {...moreProps}>{title || <i>unnamed</i>}</Typography>
        )}
      </Grid>
      {actions && (
        <Grid item>
          {actions}
        </Grid>
      )}
    </Grid>
    <Box display={'block'} flexDirection={'column'} width={'100%'}>
      {definitionName ? <LinkWithCopyToClipboard label={label} name={definitionName} itemKey={'_metainfo'} {...moreProps}/> : <LinkWithCopyToClipboard label={label} {...moreProps}/>}
      {subSectionName && <LinkWithCopyToClipboard label={'sub section'} name={subSectionName} itemKey={'_subsectionmetainfo'} {...moreProps}/>}
    </Box>
  </Compartment>
}
Title.propTypes = ({
  title: PropTypes.string,
  label: PropTypes.string,
  definitionName: PropTypes.string,
  subSectionName: PropTypes.string,
  tooltip: PropTypes.string,
  actions: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
})

export function ItemChip(props) {
  return <Chip style={{marginLeft: 8, marginBottom: 3, height: 18}} size="small" {...props} />
}
