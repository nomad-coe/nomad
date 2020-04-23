import React, { useState, useContext, useEffect } from 'react'
import PropTypes from 'prop-types'
import { makeStyles } from '@material-ui/core/styles'
import { Card, Button, Tooltip, Tabs, Tab, Paper, FormControl,
  FormGroup, Checkbox, FormControlLabel, CardContent, IconButton, FormLabel, Select, MenuItem } from '@material-ui/core'
import { useQueryParam, useQueryParams, StringParam, NumberParam } from 'use-query-params'
import SearchBar from './SearchBar'
import EntryList from './EntryList'
import DatasetList from './DatasetList'
import { DisableOnLoading } from '../api'
import { domains } from '../domains'
import PeriodicTable from './PeriodicTable'
import ReloadIcon from '@material-ui/icons/Cached'
import UploadList from './UploadsList'
import GroupList from './GroupList'
import ApiDialogButton from '../ApiDialogButton'
import SearchIcon from '@material-ui/icons/Search'
import UploadsChart from './UploadsChart'
import { Quantity } from './QuantityHistogram'
import SearchContext, { searchContext } from './SearchContext'
import {objectFilter} from '../../utils'

const resultTabs = {
  'entries': {
    label: 'Entries',
    groups: {},
    component: SearchEntryList
  },
  'groups': {
    label: 'Grouped entries',
    groups: {'dft.groups_grouped': true},
    component: SearchGroupList
  },
  'uploads': {
    label: 'Uploads',
    groups: {'uploads_grouped': true},
    component: SearchUploadList
  },
  'datasets': {
    label: 'Datasets',
    groups: {'datasets_grouped': true},
    component: SearchDatasetList
  }
}

const defaultVisalizations = {
  'elements': {
    component: ElementsVisualization,
    label: 'Elements',
    description: 'Shows data as a heatmap over the periodic table'
  },
  'users': {
    component: UsersVisualization,
    label: 'Users',
    description: 'Show statistics on user metadata'
  }
}

const useSearchStyles = makeStyles(theme => ({
  root: {
    padding: theme.spacing(3)
  }
}))

/**
 * This component shows the full search interface including result lists.
 */
export default function Search(props) {
  const {
    initialVisualizationTab,
    initialOwner,
    ownerTypes,
    initialDomain,
    initialMetric,
    initialResultTab,
    availableResultTabs,
    query,
    initialQuery,
    resultListProps,
    initialRequest,
    ...rest} = props
  const classes = useSearchStyles()
  return <DisableOnLoading>
    <SearchContext query={query} initialQuery={initialQuery}>
      <div className={classes.root} {...rest}>
        <SearchEntry
          initialTab={initialVisualizationTab}
          initialOwner={initialOwner}
          ownerTypes={ownerTypes}
          initialDomain={initialDomain}
          initialMetric={initialMetric}
          initialRequest={initialRequest}
        />
        <SearchResults
          initialTab={initialResultTab}
          availableTabs={availableResultTabs}
          resultListProps={resultListProps}
        />
      </div>
    </SearchContext>
  </DisableOnLoading>
}
Search.propTypes = {
  initialResultTab: PropTypes.string,
  initialVisualizationTab: PropTypes.string,
  availableResultTabs: PropTypes.arrayOf(PropTypes.string),
  initialOwner: PropTypes.string,
  ownerTypes: PropTypes.arrayOf(PropTypes.string),
  initialDomain: PropTypes.string,
  initialMetric: PropTypes.string,
  initialRequest: PropTypes.object,
  resultListProps: PropTypes.object,
  /**
   * Additional search parameters that will be added to all searches that are send to
   * the API. The idea is that this can be used to lock some aspects of the search for
   * special contexts, like the dataset page for example.
   */
  query: PropTypes.object,
  /**
   * Similar to query, but these parameters can be changes by the user interacting with
   * the component.
   */
  initialQuery: PropTypes.object
}

const useSearchEntryStyles = makeStyles(theme => ({
  search: {
    marginTop: theme.spacing(2),
    marginBottom: theme.spacing(8),
    maxWidth: 1024,
    margin: 'auto',
    width: '100%'
  },
  searchIcon: {
    margin: `${theme.spacing(1)}px 0`,
    padding: `6px 0 2px 0`
  },
  domainButton: {
    margin: theme.spacing(1)
  },
  metricButton: {
    margin: theme.spacing(1),
    marginRight: 0
  },
  searchBar: {
    width: '100%'
  },
  selectButton: {
    margin: theme.spacing(1)
  },
  visualizations: {
    display: 'block',
    maxWidth: 900,
    margin: 'auto',
    marginTop: theme.spacing(2),
    marginBottom: theme.spacing(2)
  }
}))
function SearchEntry({initialTab, initialOwner, ownerTypes, initialDomain, initialMetric}) {
  const classes = useSearchEntryStyles()
  const [openVisualizationParam, setOpenVisualizationParam] = useQueryParam('visualization', StringParam)
  const {domain} = useContext(searchContext)

  const visualizations = {}
  Object.assign(visualizations, defaultVisalizations)
  Object.assign(visualizations, domain.searchVisualizations)
  const openVisualizationKey = openVisualizationParam || initialTab
  const openVisualizationTab = visualizations[openVisualizationKey]

  const VisualizationComponent = openVisualizationTab ? openVisualizationTab.component : React.Fragment

  const handleVisualizationChange = value => {
    if (value === openVisualizationKey) {
      setOpenVisualizationParam('none')
    } else {
      setOpenVisualizationParam(value)
    }
  }

  return <div>
    <div className={classes.search}>
      <FormGroup row style={{alignItems: 'center'}}>
        <FormControl className={classes.searchIcon}>
          <FormLabel>
            <SearchIcon/>
          </FormLabel>
        </FormControl>
        <DomainSelect classes={{root: classes.domainButton}} initialDomain={initialDomain} />
        <OwnerSelect ownerTypes={ownerTypes} initialOwner={initialOwner}/>
        <div style={{flexGrow: 1}} />
        <VisualizationSelect
          classes={{button: classes.selectButton}}
          value={openVisualizationKey}
          onChange={handleVisualizationChange}
          visualizations={visualizations}
        />
        <MetricSelect classes={{root: classes.metricButton}} initialMetric={initialMetric}/>
      </FormGroup>

      <SearchBar classes={{autosuggestRoot: classes.searchBar}} />
    </div>

    <div className={classes.visualizations}>
      <VisualizationComponent/>
    </div>
  </div>
}
SearchEntry.propTypes = {
  initialTab: PropTypes.string,
  initialOwner: PropTypes.string,
  initialDomain: PropTypes.string,
  initialMetric: PropTypes.string,
  ownerTypes: PropTypes.arrayOf(PropTypes.string)
}

function UsersVisualization(props) {
  const {domain, response: {metric}, setStatistics} = useContext(searchContext)
  useEffect(() => {
    setStatistics(['uploader'])
  }, [])
  return <div>
    <Card>
      <CardContent>
        <UploadsChart metricsDefinitions={domain.searchMetrics}/>
      </CardContent>
    </Card>
    <Quantity quantity="uploader" title="Uploaders" scale={1} metric={metric} />
  </div>
}

function ElementsVisualization(props) {
  const [exclusive, setExclusive] = useState(false)
  const {response: {statistics, metric}, query, setQuery, setStatistics} = useContext(searchContext)
  useEffect(() => {
    setStatistics(['atoms'])
  }, [])

  const handleExclusiveChanged = () => {
    if (!exclusive) {
      setQuery({only_atoms: query.atoms, atoms: []})
    } else {
      setQuery({atoms: query.only_atoms, only_atoms: []})
    }
    setExclusive(!exclusive)
  }
  const handleAtomsChanged = atoms => {
    if (exclusive) {
      setExclusive(false)
    }
    setQuery({atoms: atoms, only_atoms: []})
  }

  return <Card>
    <CardContent>
      <PeriodicTable
        aggregations={statistics.atoms}
        metric={metric}
        exclusive={exclusive}
        values={[...(query.atoms || []), ...(query.only_atoms || [])]}
        onChanged={handleAtomsChanged}
        onExclusiveChanged={handleExclusiveChanged}
      />
    </CardContent>
  </Card>
}

const useMetricSelectStyles = makeStyles(theme => ({
  root: {
    minWidth: 130,
    paddingLeft: theme.spacing(2)
  }
}))
function MetricSelect({initialMetric}) {
  const {domain, setMetric} = useContext(searchContext)
  const [metricParam, setMetricParam] = useQueryParam('metric', StringParam)
  const metric = metricParam || initialMetric || domain.defaultSearchMetric

  useEffect(() => setMetric(metric), [metric, setMetric])

  const metricsDefinitions = domain.searchMetrics
  const classes = useMetricSelectStyles()
  const [tooltipOpen, setTooltipOpen] = useState(false)
  const handleTooltip = bool => setTooltipOpen(bool)
  return <FormControl className={classes.root}>
    <Tooltip title="Select the metric used to represent data" open={tooltipOpen}>
      <Select
        MenuProps={{
          getContentAnchorEl: null,
          anchorOrigin: {
            vertical: 'bottom',
            horizontal: 'left'
          }
        }}
        renderValue={key => {
          const metric = metricsDefinitions[key]
          return metric.shortLabel || metric.label
        }}
        value={metric}
        onChange={event => setMetricParam(event.target.value)}
        onMouseEnter={() => handleTooltip(true)}
        onMouseLeave={() => handleTooltip(false)}
        onOpen={() => handleTooltip(false)}
      >
        {Object.keys(metricsDefinitions).map(metricKey => {
          const {label, tooltip} = metricsDefinitions[metricKey]
          return (
            <MenuItem value={metricKey} key={metricKey}>
              <Tooltip title={tooltip || ''}>
                <div>{label}</div>
              </Tooltip>
            </MenuItem>
          )
        })}
      </Select>
    </Tooltip>
  </FormControl>
}
MetricSelect.propTypes = {
  initialMetric: PropTypes.string
}

function VisualizationSelect({classes, value, onChange, visualizations}) {
  return <React.Fragment>
    {Object.keys(visualizations).map(key => {
      const visualization = visualizations[key]
      return <Tooltip key={key} title={visualization.description}>
        <Button
          size="small" className={classes.button}
          color={value === key ? 'primary' : 'default'}
          onClick={() => onChange(key)}
        >
          {visualization.label}
        </Button>
      </Tooltip>
    })}
  </React.Fragment>
}
VisualizationSelect.propTypes = {
  classes: PropTypes.object.isRequired,
  value: PropTypes.string,
  visualizations: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

const useDomainSelectStyles = makeStyles(theme => ({
  root: {
    minWidth: 60,
    paddingLeft: theme.spacing(2),
    paddingRight: theme.spacing(2)
  }
}))
function DomainSelect({initialDomain}) {
  const {setDomain} = useContext(searchContext)
  const [domainParam, setDomainParam] = useQueryParam('domain', StringParam)
  const domain = domainParam || initialDomain || domains.dft.key

  useEffect(() => setDomain(domain), [domain, setDomain])

  const classes = useDomainSelectStyles()
  const [tooltipOpen, setTooltipOpen] = useState(false)
  const handleTooltip = bool => setTooltipOpen(bool)
  return <FormControl className={classes.root}>
    <Tooltip
      title="Select the data domain to search. Different domains contain different type of data."
      open={tooltipOpen}
    >
      <Select
        MenuProps={{
          getContentAnchorEl: null,
          anchorOrigin: {
            vertical: 'bottom',
            horizontal: 'left'
          }
        }}
        renderValue={key => domains[key].name}
        value={domain}
        onChange={event => setDomainParam(event.target.value)}
        onMouseEnter={() => handleTooltip(true)}
        onMouseLeave={() => handleTooltip(false)}
        onOpen={() => handleTooltip(false)}
      >
        {Object.keys(domains).map(domainKey => {
          const domain = domains[domainKey]
          return (
            <MenuItem value={domain.key} key={domain.key}>
              <Tooltip title={domain.about}>
                <div>{domain.label}</div>
              </Tooltip>
            </MenuItem>
          )
        })}
      </Select>
    </Tooltip>
  </FormControl>
}
DomainSelect.propTypes = {
  initialDomain: PropTypes.string
}

const ownerLabel = {
  all: 'All entries',
  visible: 'Include your private entries',
  public: 'Only public entries',
  user: 'Only your entries',
  staging: 'Staging area only'
}

const ownerTooltips = {
  all: 'This will show all entries in the database.',
  visible: 'Do also show entries that are only visible to you.',
  public: 'Do not show entries with embargo.',
  user: 'Do only show entries visible to you.',
  staging: 'Will only show entries that you uploaded, but not yet published.'
}

function OwnerSelect(props) {
  const {ownerTypes, initialOwner} = props
  const {setOwner} = useContext(searchContext)

  const ownerTypesToRender = ownerTypes.length === 2 ? [ownerTypes[1]] : ownerTypes

  const [ownerParam, setOwnerParam] = useQueryParam('owner', StringParam)
  const owner = ownerParam || initialOwner || 'all'

  useEffect(() => {
    setOwner(owner)
  }, [owner, setOwner])

  const handleChange = (event) => {
    if (ownerTypes.length === 2) {
      setOwnerParam(event.target.checked ? ownerTypes[1] : ownerTypes[0])
    } else {
      setOwnerParam(event.target.value)
    }
  }

  if (ownerTypes.length === 1) {
    return <React.Fragment/>
  }

  return <FormControl>
    <FormGroup row>
      {ownerTypesToRender.map(ownerToRender => (
        <Tooltip key={ownerToRender} title={ownerTooltips[ownerToRender]}>
          <FormControlLabel
            control={<Checkbox
              checked={owner === ownerToRender}
              onChange={handleChange} value="owner"
            />}
            label={ownerLabel[ownerToRender]}
          />
        </Tooltip>
      ))}
    </FormGroup>
  </FormControl>
}
OwnerSelect.propTypes = {
  ownerTypes: PropTypes.arrayOf(PropTypes.string).isRequired,
  initialOwner: PropTypes.string
}

const useSearchResultStyles = makeStyles(theme => ({
  root: theme.spacing(4)
}))
function SearchResults({availableTabs = ['entries'], initialTab = 'entries', resultListProps = {}}) {
  const classes = useSearchResultStyles()
  const {domain, setGroups} = useContext(searchContext)
  let [openTab, setOpenTab] = useQueryParam('results', StringParam)
  openTab = openTab || initialTab
  const ResultList = resultTabs[openTab].component
  const handleTabChange = tab => {
    setOpenTab(tab)
    setGroups(resultTabs[tab].groups)
  }

  useEffect(() => {
    if (openTab !== 'entries') {
      handleTabChange(openTab)
    }
  }, [])

  return <div className={classes.root}>
    <Paper>
      <Tabs
        value={openTab}
        indicatorColor="primary"
        textColor="primary"
        onChange={(event, value) => handleTabChange(value)}
      >
        {availableTabs.filter(tab => domain.searchTabs.includes(tab)).map(key => {
          const tab = resultTabs[key]
          return <Tab key={key} label={tab.label} value={key} />
        })}
      </Tabs>

      <ResultList domain={domain} {...resultListProps} />
    </Paper>
  </div>
}
SearchResults.propTypes = {
  'availableTabs': PropTypes.arrayOf(PropTypes.string),
  'initialTab': PropTypes.string,
  'resultListProps': PropTypes.object
}

function ReRunSearchButton() {
  const {update} = useContext(searchContext)
  return <Tooltip title="Re-execute the search">
    <IconButton onClick={update}>
      <ReloadIcon />
    </IconButton>
  </Tooltip>
}

const usePagination = () => {
  const {setRequestParameters} = useContext(searchContext)
  let [requestQueryParameters, setRequestQueryParameters] = useQueryParams({
    order: NumberParam, order_by: StringParam, per_page: NumberParam, page: NumberParam
  })
  requestQueryParameters = objectFilter(requestQueryParameters, key => requestQueryParameters[key])
  useEffect(
    () => setRequestParameters(requestQueryParameters),
    [requestQueryParameters, setRequestParameters]
  )
  return setRequestQueryParameters
}

const useScroll = (apiGroupName, afterParameterName) => {
  afterParameterName = afterParameterName || `${apiGroupName}_after`
  const apiAfterParameterName = `${apiGroupName}_grouped_after`

  const {response, setRequestParameters} = useContext(searchContext)
  const [queryAfterParameter, setQueryAfterParameter] = useQueryParam(afterParameterName, StringParam)
  useEffect(
    () => {
      const requestParameters = {}
      requestParameters[apiAfterParameterName] = queryAfterParameter
      setRequestParameters(requestParameters)
    }, [queryAfterParameter, setRequestParameters]
  )

  const responseGroup = response[`${apiGroupName}_grouped`]
  const after = responseGroup && responseGroup.after
  const result = {
    total: response.statistics.total.all[apiGroupName],
    onChange: requestParameters => setQueryAfterParameter(requestParameters[apiAfterParameterName])
  }
  result[afterParameterName] = after
  return result
}

function SearchEntryList(props) {
  const {response, requestParameters, apiQuery} = useContext(searchContext)
  const setRequestParameters = usePagination()
  return <EntryList
    query={apiQuery}
    editable={apiQuery.owner === 'staging' || apiQuery.owner === 'user'}
    data={response}
    onChange={setRequestParameters}
    actions={
      <React.Fragment>
        <ReRunSearchButton/>
        <ApiDialogButton data={response} />
      </React.Fragment>
    }
    {...requestParameters}
    {...props}
  />
}

function SearchDatasetList(props) {
  const {response} = useContext(searchContext)
  return <DatasetList
    data={response}
    actions={<ReRunSearchButton/>}
    {...response} {...props} {...useScroll('datasets')}
  />
}

function SearchGroupList(props) {
  const {response} = useContext(searchContext)
  return <GroupList
    data={response}
    actions={<ReRunSearchButton/>}
    {...response} {...props} {...useScroll('dft.groups', 'groups_after')}
  />
}

function SearchUploadList(props) {
  const {response} = useContext(searchContext)
  return <UploadList data={response}
    actions={<ReRunSearchButton/>}
    {...response} {...props} {...useScroll('uploads')}
  />
}
