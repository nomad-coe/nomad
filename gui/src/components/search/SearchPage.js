import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import { FormControl, FormControlLabel, Checkbox, FormGroup,
  FormLabel, IconButton, Typography, Divider, Tooltip, Tabs, Tab } from '@material-ui/core'
import { compose } from 'recompose'
import { withErrors } from '../errors'
import { withApi, DisableOnLoading } from '../api'
import SearchBar from './SearchBar'
import EntryList from './EntryList'
import SearchAggregations from './SearchAggregations'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ExpandLessIcon from '@material-ui/icons/ExpandLess'
import { withDomain } from '../domains'
import { appBase } from '../../config'
import DatasetList from './DatasetList';

export const help = `
This page allows you to **search** in NOMAD's data. The upper part of this page
gives you various options to enter and configure your search. The lower half
shows all data that fulfills your search criteria.

** Disclaimer: ** This is a preliminary version of the NOMAD software. It might
now show all of NOMAD's data. To see the full NOMAD dataset use the original
[NOMAD CoE Repository](https://repository.nomad-coe.eu/NomadRepository-1.1/search/)
for now.

#### Search Options

NOMAD's *domain-aware* search allows you to screen data by filtering based on
desired properties. This is different from basic *text-search* that traditional
search engines offer.

If you are logged-in, you can specify if you want to search among all data, publicly
available data, your own data, or just unpublished data in your [staging area](/uploads/).

The search bar allows you to specify various quantity values that you want to
see in your results. This includes *authors*, *comments*, *atom labels*, *code name*,
*system type*, *crystal system*, *basis set types*, and *XC functionals*.
Alternatively, you can click the periodic table and statistic bars to filter for respective
quantities.

The periodic table and bar-charts show metrics for all data that fit your criteria.
You can display *entries* (e.g. code runs), *unique entries*, and *datasets*.
Other more specific metrics might be available.

#### Search Results

The results table gives you a quick overview of all entries that fit your search.
You can click entries to see more details, download data, see the archive, etc.
The *raw files* tab, will show you all files that belong to the entry and offers a download
on individual, or all files. The *archive* tab, shows you the parsed data as a tree
data structure. This view is connected to NOMAD's [meta-info](${appBase}/metainfo), which acts a schema for
all parsed data. The *log* tab, will show you a log of the entry's processing.
`

class SearchPage extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    match: PropTypes.any,
    api: PropTypes.object.isRequired,
    user: PropTypes.object,
    raiseError: PropTypes.func.isRequired,
    domain: PropTypes.object,
    loading: PropTypes.number
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
    data: SearchPage.emptySearchData,
    owner: 'all',
    searchState: {
      ...SearchAggregations.defaultState
    },
    entryListState: {
      ...EntryList.defaultState
    },
    datasetListState: {
      ...DatasetList.defaultState
    },
    showDetails: true,
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
    const { owner, entryListState, datasetListState, searchState } = {...this.state, ...changes}
    const { searchValues, ...searchStateRest } = searchState
    this.setState({...changes})

    this.props.api.search({
      owner: owner,
      ...entryListState,
      ...datasetListState,
      ...searchValues,
      ...searchStateRest
    }).then(data => {
      this.setState({
        data: data || SearchPage.emptySearchData
      })
    }).catch(error => {
      if (error.name === 'NotAuthorized' && owner !== 'all') {
        this.setState({data: SearchPage.emptySearchData, owner: 'all'})
      } else {
        this.setState({data: SearchPage.emptySearchData, owner: owner})
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
    if (prevProps.api !== this.props.api) { // login/logout case, reload results
      this.update()
    } else if (prevProps.match.path !== this.props.match.path) { // navigation case
      // update if we went back to the search
      if (this.props.match.path === '/search' && this.props.match.isExact) {
        this.update()
      }
    }
  }

  handleOwnerChange(owner) {
    this.update({owner: owner})
  }

  handleClickExpand() {
    this.setState({showDetails: !this.state.showDetails})
  }

  render() {
    const { classes, user, domain, loading } = this.props
    const { data, searchState, entryListState, datasetListState, showDetails, resultTab } = this.state
    const { searchValues } = searchState
    const { pagination: { total }, statistics } = data

    const ownerLabel = {
      all: 'All entries',
      public: 'Only public entries',
      user: 'Only your entries',
      staging: 'Staging area only'
    }

    const ownerTooltips = {
      all: 'This will show all entries in the database.',
      public: 'Do not show entries that are only visible to you.',
      user: 'Do only show entries visible to you.',
      staging: 'Will only show entries that you uploaded, but not yet published.'
    }

    const withoutLogin = ['all']

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
            <div className={classes.searchEntry}>
              <FormControl>
                <FormLabel>Filter entries and show: </FormLabel>
                <FormGroup row>
                  {['all', 'public', 'user', 'staging']
                    .filter(key => user || withoutLogin.indexOf(key) !== -1)
                    .map(owner => (
                      <Tooltip key={owner} title={ownerTooltips[owner] + (user ? '' : 'You need to be logged-in for more options.')}>
                        <FormControlLabel
                          control={
                            <Checkbox checked={this.state.owner === owner} onChange={() => this.handleOwnerChange(owner)} value="owner" />
                          }
                          label={ownerLabel[owner]}
                        />
                      </Tooltip>
                    ))}
                </FormGroup>
              </FormControl>
            </div>

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

export default compose(withApi(false), withErrors, withDomain, withStyles(SearchPage.styles))(SearchPage)
