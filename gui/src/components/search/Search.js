import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import { IconButton, Divider, Tooltip, Tabs, Tab, Paper } from '@material-ui/core'
import SearchBar from './SearchBar'
import EntryList from './EntryList'
import SearchAggregations from './SearchAggregations'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ExpandLessIcon from '@material-ui/icons/ExpandLess'
import DatasetList from './DatasetList'
import SearchContext from './SearchContext'
import { DisableOnLoading } from '../api'

/**
 * Component that comprises all search views: SearchBar, SearchAggregations (aka statistics),
 * results (EntryList, DatasetList).
 */
class Search extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    query: PropTypes.object,
    showDetails: PropTypes.bool
  }

  static styles = theme => ({
    root: {
      padding: theme.spacing.unit * 3
    },
    resultsContainer: {
    },
    searchEntry: {
      minWidth: 500,
      maxWidth: 900,
      margin: 'auto',
      width: '100%'
    },
    search: {
      marginTop: theme.spacing.unit * 4,
      marginBottom: theme.spacing.unit * 8,
      display: 'flex',
      alignItems: 'center',
      minWidth: 500,
      maxWidth: 1000,
      margin: 'auto',
      width: '100%'
    },
    searchBar: {
      width: '100%'
    },
    searchDivider: {
      width: 1,
      height: 28,
      margin: theme.spacing.unit * 0.5
    },
    searchButton: {
      padding: 10
    },
    searchResults: {}
  })


  state = {
    showDetails: this.props.showDetails,
    resultTab: 'entries'
  }

  constructor(props) {
    super(props)
    this.handleClickExpand = this.handleClickExpand.bind(this)
  }

  handleClickExpand() {
    this.setState({showDetails: !this.state.showDetails})
  }

  render() {
    const { classes, query } = this.props
    const { showDetails, resultTab } = this.state

    return (
      <div className={classes.root}>
        <SearchContext query={query}>
          <DisableOnLoading>
            <div className={classes.search}>
              <SearchBar classes={{autosuggestRoot: classes.searchBar}}
                helperText="HELPER" // helperText={helperText}
              />
              <Divider className={classes.searchDivider} />
              <Tooltip title={showDetails ? 'Hide statistics' : 'Show statistics'}>
                <IconButton className={classes.searchButton} color="secondary" onClick={this.handleClickExpand}>
                  {showDetails ? <ExpandLessIcon/> : <ExpandMoreIcon/>}
                </IconButton>
              </Tooltip>
            </div>

            <div className={classes.searchEntry}>
              <SearchAggregations showDetails={showDetails}
              />
            </div>
          </DisableOnLoading>
          <Paper className={classes.resultsContainer}>
            <Tabs
              value={resultTab}
              indicatorColor="primary"
              textColor="primary"
              onChange={(event, value) => this.setState({resultTab: value})}
            >
              <Tab label="Calculations" value="entries" />
              <Tab label="Datasets" value="datasets" />
            </Tabs>

            <div className={classes.searchResults} hidden={resultTab !== 'entries'}>
              <SearchEntryList />
            </div>
            <div className={classes.searchResults} hidden={resultTab !== 'datasets'}>
              <SearchDatasetList />
            </div>
          </Paper>
        </SearchContext>
      </div>
    )
  }
}

class SearchEntryList extends React.Component {
  static contextType = SearchContext.type

  render() {
    const {state: {response, request, query}, setRequest} = this.context

    return  <EntryList
      query={query}
      editable={query.owner === 'staging' || query.owner === 'user'}
      data={response}
      onChange={setRequest}
      {...request}
    />
  }
}

class SearchDatasetList extends React.Component {
  static contextType = SearchContext.type

  render() {
    const {state: {response}, setRequest} = this.context

    return <DatasetList data={response} total={response.statistics.total.all.datasets}
      onChange={setRequest}
      {...response}
    />
  }
}

export default withStyles(Search.styles)(Search)
