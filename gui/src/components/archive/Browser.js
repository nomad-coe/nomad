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

import React, {createRef, useCallback, useContext, useEffect, useLayoutEffect, useMemo, useRef, useState} from 'react'
import PropTypes from 'prop-types'
import { Box, Button, Card, CardContent, Chip, CircularProgress, Dialog, DialogActions, DialogContent,
  DialogContentText, DialogTitle, Grid, IconButton, makeStyles, Tooltip, Typography } from '@material-ui/core'
import grey from '@material-ui/core/colors/grey'
import classNames from 'classnames'
import { useLocation, useRouteMatch, Link } from 'react-router-dom'
import { ErrorHandler } from '../ErrorHandler'
import { useApi } from '../api'
import { useErrors } from '../errors'
import NavigateIcon from '@material-ui/icons/ArrowRight'

function escapeBadPathChars(s) {
  return s.replace(/!/g, '!0').replace(/\?/g, '!1').replace(/#/g, '!2').replace(/%/g, '!3').replace(/\\/g, '!4')
}

function unescapeBadPathChars(s) {
  return s.replace(/!4/g, '\\').replace(/!3/g, '%').replace(/!2/g, '#').replace(/!1/g, '?').replace(/!0/g, '!')
}

export function formatSubSectionName(name) {
  // return name.startsWith('section_') ? name.slice(8) : name
  return name
}

/**
 * Browsers are made out of lanes. Each lane uses an adaptor that determines how to render
 * the lane contents and what adaptor is used for the next lane (depending on what is
 * selected in this lane).
 */
export class Adaptor {
  constructor(context) {
    this.context = context

    if (new.target === Adaptor) {
      throw new TypeError('Cannot construct Abstract instances directly')
    }
  }

  /**
   * If this adaptor needs to fetch some data via the API
   */
  needToFetchData() {
    return false
  }

  /**
   * A potentially asynchronous method that is called when the browser is updated if
   * the adaptor needs to fetch data
   */
  fetchData(api) {
  }

  /**
   * Called to inform adaptors that the files of an upload has been updated.
   */
  onFilesUpdated(uploadId, path) {
  }

  /**
   * A potentially asynchronous method that is used to determine the adaptor for the
   * next lane depending on the given key/url segment.
   * @returns An adaptor that is used to render the next lane.
   */
  itemAdaptor(key) {
    return null
  }

  /**
   * Renders the contents of the current lane (the lane that this adaptor represents).
   */
  render() {
    return ''
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
  }
}))
export const Browser = React.memo(function Browser({adaptor, form}) {
  const classes = useBrowserStyles()
  const rootRef = useRef()
  const outerRef = useRef()
  const innerRef = useRef()
  const { pathname } = useLocation()
  const { url } = useRouteMatch()

  const { api } = useApi()
  const { raiseError } = useErrors()

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

  const [render, setRender] = useState(0)
  const update = useCallback((lane) => {
    if (lane) {
      const index = lanes.current.indexOf(lane)
      if (index >= 0) {
        lanes.current.splice(index)
      }
    } else if (lanes.current.length > 0) {
      lanes.current.splice(-1)
    }
    // If all lanes got updated, remove all data from the root adaptor to force reload
    if (lanes.current.length === 0) {
      adaptor.data = undefined
    }
    setRender(current => current + 1)
  }, [setRender, adaptor])
  const [, setInternalRender] = useState(0)
  const internalUpdate = useCallback(() => {
    setInternalRender(current => current + 1)
  }, [setInternalRender])
  const lanes = useRef(null)

  // do no reuse the lanes if the adaptor has changed, e.g. due to updated archive data
  useEffect(() => {
    lanes.current = null
  }, [adaptor])

  useEffect(() => {
    if (!url) {
      return
    }

    const rootPath = url.endsWith('/') ? url.substring(0, url.length - 1) : url
    const segments = ['root'].concat(pathname.substring(url.length).split('/').filter(segment => segment))

    async function computeLanes() {
      const oldLanes = lanes.current
      const newLanes = []
      for (let index = 0; index < segments.length; index++) {
        const segment = unescapeBadPathChars(segments[index])
        const prev = index === 0 ? null : newLanes[index - 1]
        let lane

        if (oldLanes && oldLanes[index]?.key === segment) {
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
            fetchDataCounter: 0,
            update: update
          }
          if (!prev) {
            lane.adaptor = adaptor
          } else {
            try {
              lane.adaptor = await prev.adaptor.itemAdaptor(segment)
            } catch (error) {
              console.log(error)
              lane.error = `The item "${segment}" could not be found.`
            }
          }
        }
        if (prev) {
          prev.next = lane
        }
        if (!lane.error && lane.adaptor.needToFetchData()) {
          try {
            await lane.adaptor.fetchData(api)
            lane.fetchDataCounter += 1
          } catch (error) {
            console.log(error)
            lane.error = `Could not fetch data for "${segment}". Bad path provided?`
          }
        }
        newLanes.push(lane)
        if (lane.error) {
          break // Ignore subsequent segments/lanes
        }
      }
      lanes.current = newLanes
    }
    computeLanes().then(() => internalUpdate())
  }, [lanes, url, pathname, adaptor, render, update, internalUpdate, api, raiseError])

  const contextValue = useMemo(() => ({
    lanes: lanes,
    update: update,
    blockUntilProcessed: undefined // Will be set when creating WaitForProcessingDialog component
  }), [lanes, update])

  if (url === undefined) {
    // Can happen when navigating to another tab, possibly with the browser's back/forward buttons
    // We want to keep the cached lanes, in case the user goes back to this tab, so return immediately.
    return
  }

  return <browserContext.Provider value={contextValue}>
    {form}
    <Card>
      <CardContent>
        <div className={classes.root} ref={rootRef} >
          <div className={classes.lanesContainer} ref={outerRef} >
            <div className={classes.lanes} ref={innerRef} >
              {lanes.current && lanes.current.map((lane, index) => (
                <Lane key={`lane${index}`} lane={lane} />
              ))}
            </div>
          </div>
        </div>
        <WaitForProcessingDialog/>
      </CardContent>
    </Card>
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
    overflowY: 'scroll'
  },
  error: {
    margin: theme.spacing(1),
    minWidth: 300
  }
}))
function Lane({lane}) {
  const classes = useLaneStyles()
  const containerRef = createRef()
  const { key, adaptor, next, fetchDataCounter, error } = lane
  lane.containerRef = containerRef

  const content = useMemo(() => {
    if (error) {
      return <div className={classes.error}>
        <Typography variant="h6" color="error">ERROR</Typography>
        <Typography color="error">{error}</Typography>
      </div>
    }
    if (!adaptor) {
      return ''
    }
    return <div className={classes.root} key={`lane:${lane.path}`} data-testid={`lane${lane.index}:${lane.key}`}>
      <div className={classes.container} ref={containerRef}>
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
  }, [key, adaptor, fetchDataCounter, next?.key, classes, error])

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

export const ItemLink = React.forwardRef(function ItemLink({itemKey, ...props}, ref) {
  const lane = useContext(laneContext)
  if (lane) {
    return <Link
      {...props}
      data-testid={`item:${itemKey}`}
      to={`${lane.path}/${encodeURI(escapeBadPathChars(itemKey))}`}
    />
  } else {
    // If this is used in a non Browser context
    return ''
  }
})
ItemLink.propTypes = {
  itemKey: PropTypes.string.isRequired
}

export function ItemButton({itemKey, ...props}) {
  // ItemButtons are often used in clickable context and we should stop the propagation
  // of click events to prevent unwanted behavior.
  const handleClick = useCallback((event) => {
    event.stopPropagation()
  }, [])
  return (
    <div onClick={handleClick}>
      <IconButton {...props} component={ItemLink} itemKey={itemKey}>
        <NavigateIcon />
      </IconButton>
    </div>
  )
}
ItemButton.propTypes = {
  itemKey: PropTypes.string.isRequired
}

export function Item({children, itemKey, disabled, highlighted, icon, actions, chip}) {
  const classes = useItemStyles()
  const lane = useContext(laneContext)
  const selected = lane.next && lane.next.key
  const isSelected = itemKey && selected === itemKey
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

const useTitleStyles = makeStyles(theme => ({
  titleGridItem: {
    flexGrow: 1,
    overflow: 'hidden',
    textOverflow: 'ellipsis'
  }
}))

export function Title({title, label, tooltip, actions, ...moreProps}) {
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
  title: PropTypes.string,
  label: PropTypes.string,
  tooltip: PropTypes.string,
  actions: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
})

const statusesProcessing = ['PENDING', 'RUNNING', 'WAITING_FOR_RESULT']

function WaitForProcessingDialog() {
  const browser = useContext(browserContext)
  const { api } = useApi()
  const [jobDef, setJobDef] = useState()
  const [apiResponse, setApiResponse] = useState()
  const [apiError, setApiError] = useState()
  const [refreshing, setRefreshing] = useState(false)
  const [errorsAcknowledged, setErrorsAcknowledged] = useState(false)

  const blockUntilProcessed = useCallback(({uploadId, apiCall, apiCallText, onSuccess, onFail}) => {
    apiCall
      .then(response => {
        setApiResponse(response)
      })
      .catch(error => {
        setApiError(error)
      })
    setJobDef({uploadId, apiCallText, onSuccess, onFail})
  }, [setApiResponse, setApiError, setJobDef])

  browser.blockUntilProcessed = blockUntilProcessed // add callback to the browser context

  const refresh = useCallback(() => {
    api.get(`/uploads/${jobDef.uploadId}`)
      .then(response => {
        setApiResponse(response)
        setRefreshing(false)
      })
      .catch(error => {
        setApiError(error)
      })
  }, [api, jobDef, setApiResponse, setApiError, setRefreshing])

  const processStatus = apiResponse?.data?.process_status
  const isProcessing = jobDef && !apiError && (!apiResponse || statusesProcessing.includes(processStatus))
  const hasErrors = apiError || processStatus === 'FAILURE'
  const showDialog = jobDef && (isProcessing || (hasErrors && !errorsAcknowledged))

  useEffect(() => {
    if (jobDef) {
      if (!apiError && apiResponse && isProcessing && !refreshing) {
        // Last response was still processing. Wait and try to refresh again
        setRefreshing(true)
        const interval = setInterval(refresh(), 1000)
        return () => clearInterval(interval)
      }
      if (!showDialog) {
        // Closing dialog
        setJobDef(null)
        setApiResponse(null)
        setApiError(null)
        setRefreshing(false)
        setErrorsAcknowledged(false)
        if (hasErrors) {
          if (jobDef.onFail) {
            jobDef.onFail()
          }
        } else {
          if (jobDef.onSuccess) {
            jobDef.onSuccess()
          }
        }
      }
    }
  }, [jobDef, apiResponse, apiError, isProcessing, showDialog, hasErrors, refresh,
    refreshing, setRefreshing, setJobDef, setApiResponse, setApiError, setErrorsAcknowledged])

  if (!showDialog) {
    return ''
  }
  return <Dialog open={true} style={{textAlign: 'center'}}>
    <DialogTitle>{jobDef.apiCallText}</DialogTitle>
    <DialogContent>
      {apiResponse
        ? <DialogContentText>{apiResponse.data?.last_status_message || 'Processing...'}</DialogContentText>
        : <DialogContentText>Initiating...</DialogContentText>
      }
      {isProcessing &&
        <CircularProgress />
      }
      {apiError &&
        <DialogContentText color="error">{apiError.apiMessage || 'Operation failed'}</DialogContentText>
      }
      {processStatus === 'FAILURE' &&
        <DialogContentText color="error">{(apiResponse.data?.errors || ['Operation failed'])[0]}</DialogContentText>
      }
      {hasErrors &&
        <DialogActions>
          <Button onClick={() => setErrorsAcknowledged(true)} autoFocus>OK</Button>
        </DialogActions>
      }
    </DialogContent>
  </Dialog>
}
WaitForProcessingDialog.propTypes = {
}
