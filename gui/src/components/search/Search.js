import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import { IconButton, Typography, Divider, Tooltip, Tabs, Tab } from '@material-ui/core'
import { compose } from 'recompose'
import { withErrors } from '../errors'
import { withApi, DisableOnLoading } from '../api'
import SearchBar from './SearchBar'
import EntryList from './EntryList'
import SearchAggregations from './SearchAggregations'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ExpandLessIcon from '@material-ui/icons/ExpandLess'
import { withDomain } from '../domains'
import DatasetList from './DatasetList';
import { isEquivalent } from '../../utils';


/**
 * Component that comprises all search views: SearchBar, SearchAggregations (aka statistics),
 * results (EntryList, DatasetList).
 */
class Search extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    domain: PropTypes.object,
    loading: PropTypes.number,
    searchParameters: PropTypes.object,
    showDetails: PropTypes.bool
  }

  static styles = theme => ({
    root: {
    },
    searchContainer: {
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

  static emptySearchData = {
    results: [],
    pagination: {
      total: 0
    },
    datasets: {
      after: null,
      values: []
    },
    statistics: {
      total: {
        all: {
          datasets: 0
        }
      }
    }
  }

  state = {
    data: Search.emptySearchData,
    searchState: {
      ...SearchAggregations.defaultState
    },
    entryListState: {
      ...EntryList.defaultState
    },
    datasetListState: {
      ...DatasetList.defaultState
    },
    showDetails: this.props.showDetails,
    resultTab: 'entries'
  }

  constructor(props) {
    super(props)

    this.updateEntryList = this.updateEntryList.bind(this)
    this.updateDatasetList = this.updateDatasetList.bind(this)
    this.updateSearch = this.updateSearch.bind(this)
    this.handleClickExpand = this.handleClickExpand.bind(this)

    this._mounted = false
  }

  updateEntryList(changes) {
    const entryListState = {
      ...this.state.entryListState, ...changes
    }
    this.update({entryListState: entryListState})
  }

  updateDatasetList(changes) {
    const datasetListState = {
      ...this.state.datasetListState, ...changes
    }
    this.update({datasetListState: datasetListState})
  }

  updateSearch(changes) {
    const searchState = {
      ...this.state.searchState, ...changes
    }
    this.update({searchState: searchState})
  }

  update(changes) {
    if (!this._mounted) {
      return
    }

    changes = changes || {}
    const { searchParameters } = this.props
    const { entryListState, datasetListState, searchState } = {...this.state, ...changes}
    const { searchValues, ...searchStateRest } = searchState
    this.setState({...changes})

    this.props.api.search({
      datasets: true,
      statistics: true,
      ...entryListState,
      ...datasetListState,
      ...searchValues,
      ...searchStateRest,
      ...searchParameters
    }).then(data => {
      this.setState({
        data: data || Search.emptySearchData
      })
    }).catch(error => {
      if (error.name === 'NotAuthorized' && this.props.searchParameters.owner !== 'all') {
        this.setState({data: Search.emptySearchData})
      } else {
        this.setState({data: Search.emptySearchData})
        this.props.raiseError(error)
      }
    })
  }

  componentDidMount() {
    this._mounted = true
    this.update()
  }

  componentWillUnmount() {
    this._mounted = false
  }

  componentDidUpdate(prevProps) {
    // login/logout or changed search paraemters -> reload results
    if (prevProps.api !== this.props.api || !isEquivalent(prevProps.searchParameters, this.props.searchParameters)) {
      this.update()
    }
  }

  handleClickExpand() {
    this.setState({showDetails: !this.state.showDetails})
  }

  render() {
    const { classes, domain, loading } = this.props
    const { data, searchState, entryListState, datasetListState, showDetails, resultTab } = this.state
    const { searchValues } = searchState
    const { pagination: { total }, statistics } = data

    const helperText = <span>
      There are {Object.keys(domain.searchMetrics).filter(key => statistics.total.all[key]).map(key => {
        return <span key={key}>
          {domain.searchMetrics[key].renderResultString(!loading && statistics.total.all[key] !== undefined ? statistics.total.all[key] : '...')}
        </span>
      })}{Object.keys(searchValues).length ? ' left' : ''}.
    </span>

    return (
      <div className={classes.root}>
        <div className={classes.searchContainer}>
          <DisableOnLoading>
            <div className={classes.search}>
              <SearchBar classes={{autosuggestRoot: classes.searchBar}}
                fullWidth fullWidthInput={false} helperText={helperText}
                label="search"
                placeholder={domain.searchPlaceholder}
                data={data} searchValues={searchValues}
                InputLabelProps={{
                  shrink: true
                }}
                onChanged={values => this.updateSearch({searchValues: values})}
              />
              <Divider className={classes.searchDivider} />
              <Tooltip title={showDetails ? 'Hide statistics' : 'Show statistics'}>
                <IconButton className={classes.searchButton} color="secondary" onClick={this.handleClickExpand}>
                  {showDetails ? <ExpandLessIcon/> : <ExpandMoreIcon/>}
                </IconButton>
              </Tooltip>
            </div>

            <div className={classes.searchEntry}>
              <SearchAggregations
                data={data} {...searchState} onChange={this.updateSearch}
                showDetails={showDetails}
              />
            </div>
          </DisableOnLoading>
        </div>
        <div className={classes.resultsContainer}>
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
            <Typography variant="caption" style={{margin: 12}}>
              About {total.toLocaleString()} results:
            </Typography>

            <EntryList
              data={data} total={total}
              onChange={this.updateEntryList}
              {...entryListState}
            />
          </div>
          <div className={classes.searchResults} hidden={resultTab !== 'datasets'}>
            <Typography variant="caption" style={{margin: 12}}>
              About {statistics.total.all.datasets.toLocaleString()} datasets:
            </Typography>

            <DatasetList data={data} total={statistics.total.all.datasets}
              onChange={this.updateDatasetList}
              {...datasetListState}
            />
          </div>
        </div>
      </div>
    )
  }
}

export default compose(withApi(false), withErrors, withDomain, withStyles(Search.styles))(Search)
