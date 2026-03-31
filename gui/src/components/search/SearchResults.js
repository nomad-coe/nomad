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
import React, { useState, useMemo, useEffect, useCallback } from 'react'
import PropTypes from 'prop-types'
import jmespath from 'jmespath'
import {Paper, Typography, Tooltip, IconButton, Button, Dialog, DialogTitle, DialogContent} from '@material-ui/core'
import { Alert } from '@material-ui/lab'
import Icon from '@material-ui/core/Icon'
import GitHubIcon from '@material-ui/icons/GitHub'
import {
  Datatable,
  DatatableLoadMorePagination,
  DatatableTable,
  DatatableToolbar,
  DatatableToolbarActions
} from '../datatable/Datatable'
import EntryDownloadButton from '../entry/EntryDownloadButton'
import EntryDetails, { EntryRowActions } from '../entry/EntryDetails'
import { MaterialRowActions } from '../material/MaterialDetails'
import { pluralize, formatInteger, parseJMESPath } from '../../utils'
import { isEmpty, isArray } from 'lodash'
import { useSearchContext } from './SearchContext'
import {useTools} from "../north/NorthPage"
import DialogActions from "@material-ui/core/DialogActions"
import {getIconUrl, useNorthToolHook} from "../north/NorthTool"
import {useKeycloak} from "@react-keycloak/web"
import {useErrors} from "../errors"
import {northBase} from "../../config"

/**
 * Helper function to resolve a jmespath-like path to a quantity in an entry
 */
export const extractPathContent = ({action, data, key}) => {
  const {path} = parseJMESPath(action[key])
  return jmespath.search(data, path)
}

/**
 * Used to retrieve an URL link from the row metadata and display a link icon to
 * that resource.
 */
export const ActionURL = React.memo(({action, data}) => {
  const {raiseError} = useErrors()
  let href = extractPathContent({action, data, key: 'path'})

  if (isArray(href)) {
    if (href.length > 1) {
      raiseError(`
        The data for a row action as returned by the path "${action.path}"
        contains an array with more than one item. Update the path in your app
        to target only one item in the array.`
      )
      return null
    } else {
      href = href[0]
    }
  }
  const disabled = !href
  const size = 'medium'
  const svgIcon = {
    'github': <GitHubIcon fontSize={size}/>
  }[action.icon]

  return <Tooltip title={disabled ? 'Not available' : (action.description || '')}>
    <div>
      <IconButton size={size} href={href} target="_blank" disabled={disabled}>
        {svgIcon || <Icon fontSize={size}>{action?.icon || 'launch'}</Icon>}
      </IconButton>
    </div>
  </Tooltip>
})
ActionURL.propTypes = {
  action: PropTypes.object.isRequired, // Action configuration from app config
  data: PropTypes.object.isRequired // ES index data
}

/**
 * Used to build a north container and open it in a new tab.
 */
export const NorthURL = React.memo(({ action, data }) => {
  const {raiseError} = useErrors()
  const { keycloak } = useKeycloak()
  const isAuthenticated = keycloak.authenticated

  const northTools = useTools()
  const [openDialog, setOpenDialog] = useState(false)
  const [toolUrl, setToolUrl] = useState(undefined)
  const [dialogMessage, setDialogMessage] = useState("")
  const [showStopOnly, setShowStopOnly] = useState(false)

  const filepath = extractPathContent({ action, data, key: "filepath" })

  const toolsData = action?.tool_name
    ? Object.keys(northTools)
        .filter(tool => tool === action.tool_name)
        .map(key => {
          const { with_path, mount_path, icon, path_prefix } = northTools[key].north_tool
          return {
            name: key,
            title: key,
            with_path,
            mount_path,
            icon,
            path_prefix
          }
        })
    : []

  const tool = toolsData.length > 0 ? toolsData[0] : null
  const { launch, stop, getToolStatus } = useNorthToolHook(tool ?? {}, data?.upload_id ?? "", filepath)

  const handleLaunch = async () => {
    if (!tool) return
    const { state: currentState, uploadid_is_mounted } = await getToolStatus()
    if (currentState === "stopped") {
      const { toolUrl: newUrl } = await launch()
      setToolUrl(newUrl)
      if (newUrl) window.open(newUrl, tool.name)
    } else if (uploadid_is_mounted === false) {
      setDialogMessage(
        `The upload ${data?.upload_id} is not mounted. You need to stop the tool before relaunching it.
        Stopping the tool may cause loss of unsaved data.`)
      setShowStopOnly(true)
      setOpenDialog(true)
    } else {
      setDialogMessage("The tool is already running. Stopping it will delete any unsaved data.")
      setShowStopOnly(false)
      setOpenDialog(true)
    }
  }

  const handleStopAndRetry = async () => {
    await handleStop()
    setOpenDialog(false)
  }

  const handleStop = async () => {
    if (!tool) return
    await stop()
    setOpenDialog(false)
  }

  const handleOpenTool = async () => {
  setOpenDialog(false)
  let url = toolUrl

  if (!url) {
    const { state, uploadid_is_mounted } = await getToolStatus()
    if (state === "running" && uploadid_is_mounted) {
      const { toolUrl: fetchedUrl } = await launch()
      url = fetchedUrl || `${northBase}/user/${keycloak.tokenParsed?.preferred_username}/${tool.name}`
    } else {
      url = `${northBase}/user/${keycloak.tokenParsed?.preferred_username}/${tool.name}`
    }
  }

  window.open(url, "_blank")
}

  const disabled = !filepath || !isAuthenticated
  const tooltipMessage = !isAuthenticated
    ? "You must be logged in to use this tool."
    : !filepath
    ? "No file path found. Tool cannot be launched."
    : action.description || ""

  if (!tool) {
    return null
  }

  if (isArray(filepath)) {
    raiseError(`
    Encountered and array in ${action.filepath}. Expected a string.
    Update the path in you app to target only one item in the array.`
    )
    return
  }

  return (
    <>
      <Tooltip title={tooltipMessage}>
        <div>
          <IconButton
            onClick={handleLaunch}
            size="medium"
            disabled={disabled}
            style={{ filter: disabled ? "grayscale(100%)" : "none" }}
          >
            <Icon fontSize="medium">
              {tool.icon ? (
                <img
                  style={{ width: "30px", height: "30px", objectFit: "contain" }}
                  src={getIconUrl(tool.icon)}
                  alt="icon"
                />
              ) : (
                <Icon fontSize="medium">{action?.icon || "launch"}</Icon>
              )}
            </Icon>
          </IconButton>
        </div>
      </Tooltip>

      <Dialog open={openDialog} onClose={() => setOpenDialog(false)}>
        <DialogTitle>Tool Status</DialogTitle>
        <DialogContent>
          <p>{dialogMessage}</p>
        </DialogContent>
        <DialogActions>
          <Button onClick={handleStopAndRetry} color="error">
            Stop {action.tool_name}
          </Button>
          {!showStopOnly && (
            <Button onClick={handleOpenTool} color="primary">
              Open {action.tool_name}
            </Button>
          )}
          <Button onClick={() => setOpenDialog(false)}>Cancel</Button>
        </DialogActions>
      </Dialog>
    </>
  )
})
NorthURL.propTypes = {
  action: PropTypes.object.isRequired,
  data: PropTypes.object.isRequired
}

/**
 * Displays the list of search results.
 */
export const SearchResults = React.memo(({
  noAction,
  onSelectedChanged,
  defaultUncollapsedEntryID,
  title,
  'data-testid': testID,
  PaperProps,
  ...otherProps
}) => {
  const {columns, resource, rows, useResults, useApiQuery} = useSearchContext()
  const {data, pagination, setPagination} = useResults()
  const apiQuery = useApiQuery()
  const [selected, setSelected] = useState(new Set())
  const shownColumns = columns
    ? columns
      .filter((column) => column.selected)
      .map((column) => column.search_quantity)
    : []

  useEffect(() => {
    if (onSelectedChanged) {
      onSelectedChanged(selected)
    }
  }, [onSelectedChanged, selected])

  const query = useMemo(() => {
    if (selected === 'all') {
      return apiQuery
    }
    return {'entry_id:any': [...selected]}
  }, [selected, apiQuery])

  const actions = useCallback((data) => {
    const actionComponents = [];
    (rows?.actions?.items || []).forEach((action, i) => {
      const component = {
        url: <ActionURL key={i} action={action} data={data}/>,
        north: <NorthURL key={i} action={action} data={data}/>
      }[action.type]
      if (component) {
        actionComponents.push(component)
      }
    })
    if (resource === "entries") {
      actionComponents.push(<EntryRowActions key="entries-default-action" data={data}/>)
    } else if (resource === "materials") {
      actionComponents.push(<MaterialRowActions key="materials-default-action" data={data}/>)
    }
    return actionComponents
  }, [rows, resource])

  if (isEmpty(columns)) {
    return <Alert severity="warning">
      No search columns defined within this search context. Ensure that all GUI artifacts are created.
    </Alert>
  }

  if (pagination.total === 0) {
    return <Typography>no results</Typography>
  }

  if (!pagination.total) {
    return <Typography>searching ...</Typography>
  }

  // Select components based on the targeted resource
  let details
  let buttons
  if (resource === "entries") {
    details = rows?.details?.render || EntryDetails
    if (!noAction) buttons = <EntryDownloadButton tooltip="Download files" query={query} />
  }

  return <Paper data-testid={testID} {...PaperProps}>
    <Datatable
      data={data}
      pagination={pagination}
      onPaginationChanged={setPagination}
      columns={columns}
      shownColumns={shownColumns}
      selected={rows?.selection?.enabled ? selected : undefined}
      getId={option => option.entry_id}
      onSelectedChanged={rows?.selection?.enabled ? setSelected : undefined}
      {...otherProps}
    >
      <DatatableToolbar title={`${formatInteger(data.length)}/${pluralize(title || 'result', pagination.total, true, true, title ? undefined : 'search')}`}>
        {rows?.selection?.enabled &&
          <DatatableToolbarActions selection>
            {buttons}
          </DatatableToolbarActions>
        }
      </DatatableToolbar>
      <DatatableTable
        actions={rows?.actions?.enabled ? actions : undefined}
        details={rows?.details?.enabled ? details : undefined}
        defaultUncollapsedRow={defaultUncollapsedEntryID && data.find(row => row.entry_id === defaultUncollapsedEntryID)}
      >
        <DatatableLoadMorePagination color="primary">load more</DatatableLoadMorePagination>
      </DatatableTable>
    </Datatable>
  </Paper>
})

SearchResults.propTypes = {
  noAction: PropTypes.bool,
  PaperProps: PropTypes.object,
  onSelectedChanged: PropTypes.func,
  defaultUncollapsedEntryID: PropTypes.string,
  title: PropTypes.string,
  'data-testid': PropTypes.string
}

SearchResults.defaultProps = {
  'data-testid': 'search-results'
}
