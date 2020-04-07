import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import { Card, Button, List, ListItem, ListItemText, Tooltip, Tabs, Tab, Paper, FormControl,
  FormGroup, Checkbox, FormControlLabel, Popover, CardContent, IconButton, FormLabel } from '@material-ui/core'
import SearchBar from './SearchBar'
import EntryList from './EntryList'
import DatasetList from './DatasetList'
import SearchContext from './SearchContext'
import { DisableOnLoading } from '../api'
import { domains } from '../domains'
import KeepState from '../KeepState'
import PeriodicTable from './PeriodicTable'
import ReloadIcon from '@material-ui/icons/Cached'
import UploadList from './UploadsList'
import GroupList from './GroupList'
import ApiDialogButton from '../ApiDialogButton'
import SearchIcon from '@material-ui/icons/Search'
import UploadsChart from './UploadsChart'
import UploadersList from './UploadersList'

class Search extends React.Component {
  static tabs = {
    'entries': {
      label: 'Entries',
      render: (props) => <SearchEntryList {...(props || {})}/>
    },
    'groups': {
      label: 'Grouped entries',
      render: (props) => <SearchGroupList {...props} />
    },
    'uploads': {
      label: 'Uploads',
      render: (props) => <SearchUploadList {...props}/>
    },
    'datasets': {
      label: 'Datasets',
      render: (props) => <SearchDatasetList {...props}/>
    }
  }

  static propTypes = {
    classes: PropTypes.object.isRequired,
    resultTab: PropTypes.string,
    entryListProps: PropTypes.object,
    visualization: PropTypes.string,
    tabs: PropTypes.arrayOf(PropTypes.string)
  }

  static styles = theme => ({
    root: {
      padding: theme.spacing.unit * 3
    },
    search: {
      marginTop: theme.spacing.unit * 2,
      marginBottom: theme.spacing.unit * 8,
      maxWidth: 1024,
      margin: 'auto',
      width: '100%'
    },
    searchIcon: {
      margin: `${theme.spacing.unit}px 0`,
      padding: `6px 0 2px 0`
    },
    domainButton: {
      margin: theme.spacing.unit
    },
    metricButton: {
      margin: theme.spacing.unit,
      marginRight: 0
    },
    searchBar: {
      width: '100%'
    },
    selectButton: {
      margin: theme.spacing.unit
    },
    visalizations: {
      display: 'block',
      maxWidth: 900,
      margin: 'auto',
      marginTop: theme.spacing.unit * 2,
      marginBottom: theme.spacing.unit * 2
    },
    searchResults: {
      marginTop: theme.spacing.unit * 4
    }
  })

  static visalizations = {
    'elements': {
      render: props => <ElementsVisualization {...props}/>,
      label: 'Elements',
      description: 'Shows data as a heatmap over the periodic table'
    },
    'domain': {
      render: props => <DomainVisualization {...props}/>,
      label: 'Meta data',
      description: 'Shows histograms on key metadata'
    },
    'property': {
      render: props => <PropertyVisualization {...props}/>,
      label: 'Properties',
      description: 'Shows histograms on key properties'
    },
    'users': {
      render: props => <UsersVisualization {...props}/>,
      label: 'Users',
      description: 'Show statistics on user metadata'
    }
  }

  static contextType = SearchContext.type

  state = {
    resultTab: 'entries',
    openVisualization: this.props.visualization
  }

  constructor(props) {
    super(props)
    this.handleVisualizationChange = this.handleVisualizationChange.bind(this)
  }

  handleVisualizationChange(value) {
    const {openVisualization} = this.state
    if (value === openVisualization) {
      this.setState({openVisualization: null})
    } else {
      this.setState({openVisualization: value})
    }
  }

  componentDidMount() {
    if ((this.props.resultTab || 'entries') !== 'entries') {
      this.setState({resultTab: this.props.resultTab})
    }
  }

  handleTabChange(tab) {
    const {setRequest} = this.context

    this.setState({resultTab: tab}, () => {
      setRequest({
        uploads_grouped: tab === 'uploads' ? true : undefined,
        datasets_grouped: tab === 'datasets' ? true : undefined,
        'dft.groups_grouped': tab === 'groups' ? true : undefined
      })
    })
  }

  render() {
    const {classes, entryListProps, tabs} = this.props
    const {resultTab, openVisualization} = this.state
    const {domain} = this.context.state
    // const {state: {request: {uploads, datasets, groups}}} = this.context

    return <DisableOnLoading>
      <div className={classes.root}>
        <div className={classes.search}>
          <FormGroup row>
            <FormControl className={classes.searchIcon}><FormLabel><SearchIcon/></FormLabel></FormControl>
            <DomainSelect classes={{root: classes.domainButton}} />
            <OwnerSelect />
            <div style={{flexGrow: 1}} />
            <VisualizationSelect
              classes={{button: classes.selectButton}}
              value={openVisualization}
              onChange={this.handleVisualizationChange}
            />
            <MetricSelect classes={{root: classes.metricButton}} />
          </FormGroup>

          <SearchBar classes={{autosuggestRoot: classes.searchBar}} />
        </div>

        <div className={classes.visalizations}>
          {Object.keys(Search.visalizations).map(key => {
            return Search.visalizations[key].render({
              key: key, open: openVisualization === key
            })
          })
          }
        </div>

        <div className={classes.searchResults}>
          <Paper>
            <Tabs
              value={resultTab}
              indicatorColor="primary"
              textColor="primary"
              onChange={(event, value) => this.handleTabChange(value)}
            >
              {tabs.filter(tab => domain.searchTabs.includes(tab)).map(tab => <Tab
                key={tab}
                label={Search.tabs[tab].label}
                value={tab}
              />)}
            </Tabs>

            {tabs.map(tab => <KeepState
              key={tab}
              visible={resultTab === tab}
              render={() => Search.tabs[tab].render({domain: domain, ...entryListProps})}
            />)}
          </Paper>
        </div>
      </div>
    </DisableOnLoading>
  }
}

class DomainVisualization extends React.Component {
  static propTypes = {
    open: PropTypes.bool
  }

  static contextType = SearchContext.type

  render() {
    const {domain} = this.context.state
    const {open} = this.props

    return <KeepState visible={open} render={() =>
      <domain.SearchAggregations />
    }/>
  }
}

class PropertyVisualization extends React.Component {
  static propTypes = {
    open: PropTypes.bool
  }

  static contextType = SearchContext.type

  render() {
    const {domain} = this.context.state
    const {open} = this.props

    return <KeepState visible={open} render={() =>
      <domain.SearchByPropertyAggregations />
    }/>
  }
}

class UsersVisualization extends React.Component {
  static propTypes = {
    open: PropTypes.bool
  }

  static contextType = SearchContext.type

  render() {
    const {domain} = this.context.state
    const {open} = this.props

    return <KeepState visible={open} render={() =>
      <div>
        <Card>
          <CardContent>
            <UploadsChart metricsDefinitions={domain.searchMetrics}/>
          </CardContent>
        </Card>
        <UploadersList />
      </div>
    }/>
  }
}

class ElementsVisualization extends React.Component {
  static propTypes = {
    open: PropTypes.bool
  }

  static contextType = SearchContext.type

  constructor(props) {
    super(props)
    this.handleExclusiveChanged = this.handleExclusiveChanged.bind(this)
    this.handleAtomsChanged = this.handleAtomsChanged.bind(this)
  }

  state = {
    exclusive: false
  }

  handleExclusiveChanged() {
    this.setState({exclusive: !this.state.exclusive}, () => {
      const {state: {query}, setQuery} = this.context
      if (this.state.exclusive) {
        setQuery({...query, only_atoms: query['atoms'], atoms: []})
      } else {
        setQuery({...query, atoms: query.only_atoms, only_atoms: []})
      }
    })
  }

  handleAtomsChanged(atoms) {
    if (this.state.exclusive) {
      this.setState({exclusive: false})
    }

    const {state: {query}, setQuery} = this.context
    setQuery({...query, atoms: atoms, only_atoms: []})
  }

  render() {
    const {open} = this.props
    const {state: {response: {statistics}, query, metric}} = this.context

    return <KeepState visible={open} render={() =>
      <Card>
        <CardContent>
          <PeriodicTable
            aggregations={statistics.atoms}
            metric={metric}
            exclusive={this.state.exclusive}
            values={[...(query.atoms || []), ...(query.only_atoms || [])]}
            onChanged={this.handleAtomsChanged}
            onExclusiveChanged={this.handleExclusiveChanged}
          />
        </CardContent>
      </Card>
    }/>
  }
}

class MetricSelect extends React.Component {
  static contextType = SearchContext.type

  constructor(props) {
    super(props)
    this.handleClick = this.handleClick.bind(this)
    this.handleClose = this.handleClose.bind(this)
  }

  state = {
    anchorEl: null
  }

  handleClick = event => {
    this.setState({
      anchorEl: event.currentTarget
    })
  }

  handleClose = () => {
    this.setState({
      anchorEl: null
    })
  }

  handleToggle = (value) => {
    const {setMetric} = this.context
    this.setState({anchorEl: null})
    setMetric(value)
  }

  render() {
    const {metric, domain} = this.context.state
    const {anchorEl} = this.state

    const metricsDefinitions = domain.searchMetrics
    const {label, shortLabel} = metricsDefinitions[metric]
    return <React.Fragment>
      <Tooltip title="Select the metric used to represent data">
        <Button size="small" onClick={this.handleClick} {...this.props} >
          {shortLabel || label} &#9662;
        </Button>
      </Tooltip>
      <Popover
        open={anchorEl !== null}
        anchorEl={anchorEl}
        onClose={this.handleClose}
        anchorOrigin={{
          vertical: 'bottom',
          horizontal: 'center'
        }}
        transformOrigin={{
          vertical: 'top',
          horizontal: 'center'
        }}
      >
        <List>
          {Object.keys(metricsDefinitions).map(key => {
            const {label, tooltip} = metricsDefinitions[key]
            return (
              <ListItem
                key={key} role={undefined} dense button
                onClick={() => this.handleToggle(key)}
              >
                <Tooltip title={tooltip || ''}>
                  <ListItemText primary={label} />
                </Tooltip>
              </ListItem>
            )
          })}
        </List>
      </Popover>
    </React.Fragment>
  }
}

class VisualizationSelect extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    value: PropTypes.string,
    onChange: PropTypes.func.isRequired
  }

  render() {
    const {classes, value, onChange} = this.props
    return <React.Fragment>
      {Object.keys(Search.visalizations).map(key => {
        const visualization = Search.visalizations[key]
        return <Tooltip key={key} title={visualization.description}>
          <Button
            size="small" variant="outlined" className={classes.button}
            color={value === key ? 'primary' : null}
            onClick={() => onChange(key)}
          >
            {visualization.label}
          </Button>
        </Tooltip>
      })}
    </React.Fragment>
  }
}

class DomainSelect extends React.Component {
  static contextType = SearchContext.type

  constructor(props) {
    super(props)
    this.handleClick = this.handleClick.bind(this)
    this.handleClose = this.handleClose.bind(this)
  }

  state = {
    anchorEl: null
  }

  handleClick = event => {
    this.setState({
      anchorEl: event.currentTarget
    })
  }

  handleClose = () => {
    this.setState({
      anchorEl: null
    })
  }

  handleToggle = (value) => {
    const { setDomain } = this.context
    setDomain(value)
    this.handleClose()
  }

  render() {
    const {domain} = this.context.state
    const {anchorEl} = this.state

    return <React.Fragment>
      <Tooltip title="Select the domain to search in">
        <Button size="small" onClick={this.handleClick} {...this.props}>
          {domain.name} &#9662;
        </Button>
      </Tooltip>
      <Popover
        open={anchorEl !== null}
        anchorEl={anchorEl}
        onClose={this.handleClose}
        anchorOrigin={{
          vertical: 'bottom',
          horizontal: 'center'
        }}
        transformOrigin={{
          vertical: 'top',
          horizontal: 'center'
        }}
      >
        <List>
          {Object.keys(domains).map(key => {
            const {label, searchTooltip} = domains[key]
            return (
              <ListItem
                key={key} role={undefined} dense button
                onClick={() => this.handleToggle(key)}
              >
                <Tooltip title={searchTooltip || ''}>
                  <ListItemText primary={label} />
                </Tooltip>
              </ListItem>
            )
          })}
        </List>
      </Popover>
    </React.Fragment>
  }
}

class OwnerSelect extends React.Component {
  static ownerLabel = {
    all: 'All entries',
    visible: 'Include your private entries',
    public: 'Only public entries',
    user: 'Only your entries',
    staging: 'Staging area only'
  }

  static ownerTooltips = {
    all: 'This will show all entries in the database.',
    visible: 'Do also show entries that are only visible to you.',
    public: 'Do not entries with embargo.',
    user: 'Do only show entries visible to you.',
    staging: 'Will only show entries that you uploaded, but not yet published.'
  }

  static contextType = SearchContext.type

  constructor(props) {
    super(props)
    this.handleChange = this.handleChange.bind(this)
  }

  handleChange(event) {
    const {props: {ownerTypes}, setQuery} = this.context
    if (ownerTypes.length === 2) {
      setQuery({owner: event.target.checked ? ownerTypes[1] : ownerTypes[0]})
    } else {
      setQuery({owner: event.target.value})
    }
  }

  render() {
    const {state: {query: {owner}}, props: {ownerTypes}} = this.context
    const selectedOwner = owner

    if (ownerTypes.length === 1) {
      return <React.Fragment/>
    }

    const ownerTypesToRender = ownerTypes.length === 2 ? [ownerTypes[1]] : ownerTypes

    return <FormControl>
      <FormGroup row>
        {ownerTypesToRender.map(owner => (
          <Tooltip key={owner} title={OwnerSelect.ownerTooltips[owner]}>
            <FormControlLabel
              control={<Checkbox
                checked={selectedOwner === owner}
                onChange={this.handleChange} value="owner"
              />}
              label={OwnerSelect.ownerLabel[owner]}
            />
          </Tooltip>
        ))}
      </FormGroup>
    </FormControl>
  }
}

class ReRunSearchButton extends React.PureComponent {
  static contextType = SearchContext.type

  render() {
    const {setRequest} = this.context

    return <Tooltip title="Re-execute the search">
      <IconButton onClick={() => setRequest({})}>
        <ReloadIcon />
      </IconButton>
    </Tooltip>
  }
}

class SearchEntryList extends React.Component {
  static contextType = SearchContext.type

  render() {
    const {state: {response, request, query}, props, setRequest} = this.context

    return <EntryList
      query={{...query, ...props.query}}
      editable={query.owner === 'staging' || query.owner === 'user'}
      data={response}
      onChange={setRequest}
      actions={
        <React.Fragment>
          <ReRunSearchButton/>
          <ApiDialogButton data={response} />
        </React.Fragment>
      }
      {...request}
      {...this.props}
    />
  }
}

class SearchDatasetList extends React.Component {
  static contextType = SearchContext.type

  render() {
    const {state: {response}, setRequest} = this.context

    return <DatasetList data={response}
      total={response.statistics.total.all.datasets}
      datasets_after={response.datasets_grouped && response.datasets_grouped.after}
      onChange={setRequest}
      actions={<ReRunSearchButton/>}
      {...response} {...this.props}
    />
  }
}

class SearchGroupList extends React.Component {
  static contextType = SearchContext.type

  render() {
    const {state: {response}, setRequest} = this.context

    return <GroupList data={response}
      total={response.statistics.total.all['dft.groups']}
      groups_after={response['dft.groups_grouped'] && response['dft.groups_grouped'].after}
      onChange={setRequest}
      actions={<ReRunSearchButton/>}
      {...response} {...this.props}
    />
  }
}

class SearchUploadList extends React.Component {
  static contextType = SearchContext.type

  render() {
    const {state: {response}, setRequest} = this.context

    return <UploadList data={response}
      total={response.statistics.total.all.uploads}
      uploads_after={response.uploads_grouped && response.uploads_grouped.after}
      onChange={setRequest}
      actions={<ReRunSearchButton/>}
      {...response} {...this.props}
    />
  }
}

export default withStyles(Search.styles)(Search)
