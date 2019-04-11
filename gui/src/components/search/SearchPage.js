import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import { FormControl, FormControlLabel, Checkbox, FormGroup,
  FormLabel, IconButton, Typography, Divider, Tooltip } from '@material-ui/core'
import { compose } from 'recompose'
import { withErrors } from '../errors'
import { withApi, DisableOnLoading } from '../api'
import SearchBar from './SearchBar'
import SearchResultList from './SearchResultList'
import SearchAggregations from './SearchAggregations'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ExpandLessIcon from '@material-ui/icons/ExpandLess'
import { withDomain } from '../domains'

class SearchPage extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    user: PropTypes.object,
    raiseError: PropTypes.func.isRequired,
    domain: PropTypes.object,
    loading: PropTypes.number
  }

  static styles = theme => ({
    root: {
      padding: theme.spacing.unit * 3
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
    data: {
      results: [],
      pagination: {
        total: 0
      },
      aggregations: {},
      metrics: {}
    },
    owner: 'all',
    searchState: {
      ...SearchAggregations.defaultState
    },
    searchResultListState: {
      ...SearchResultList.defaultState
    },
    showDetails: true
  }

  constructor(props) {
    super(props)

    this.updateSearchResultList = this.updateSearchResultList.bind(this)
    this.updateSearch = this.updateSearch.bind(this)
    this.handleClickExpand = this.handleClickExpand.bind(this)
  }

  updateSearchResultList(changes) {
    const searchResultListState = {
      ...this.state.searchResultListState, ...changes
    }
    this.update({searchResultListState: searchResultListState})
  }

  updateSearch(changes) {
    const searchState = {
      ...this.state.searchState, ...changes
    }
    this.update({searchState: searchState})
  }

  update(changes) {
    changes = changes || {}
    const { owner, searchResultListState, searchState } = {...this.state, ...changes}
    const { searchValues, ...searchStateRest } = searchState
    this.setState({...changes})

    this.props.api.search({
      owner: owner,
      ...searchResultListState,
      ...searchValues,
      ...searchStateRest
    }).then(data => {
      this.setState({
        data: data
      })
    }).catch(errors => {
      this.setState({data: [], total: 0, owner: owner})
      this.props.raiseError(errors)
    })
  }

  componentDidMount() {
    this.update()
  }

  componentDidUpdate(prevProps) {
    if (prevProps.api !== this.props.api) {
      this.update()
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
    const { data, searchState, searchResultListState, showDetails } = this.state
    const { searchValues } = searchState
    const { pagination: { total }, metrics } = data

    const ownerLabel = {
      all: 'All entries',
      public: 'Only public entries',
      user: 'Only your entries',
      staging: 'Only entries from your staging area'
    }

    const useMetric = Object.keys(metrics).find(metric => metric !== 'code_runs') || 'code_runs'
    const helperText = <span>
      There are {Object.keys(domain.searchMetrics).map(key => {
        return (key === useMetric || key === 'code_runs') ? <span key={key}>
          {domain.searchMetrics[key].renderResultString(!loading && metrics[key] ? metrics[key] : '...')}
        </span> : ''
      })}{Object.keys(searchValues).length ? ' left' : ''}.
    </span>

    return (
      <div className={classes.root}>
        <DisableOnLoading>
          { user
            ? <div className={classes.searchEntry}>
              <FormControl>
                <FormLabel>Filter entries and show: </FormLabel>
                <FormGroup row>
                  {['all', 'public', 'user', 'staging'].map(owner => (
                    <FormControlLabel key={owner}
                      control={
                        <Checkbox checked={this.state.owner === owner} onChange={() => this.handleOwnerChange(owner)} value="owner" />
                      }
                      label={ownerLabel[owner]}
                    />
                  ))}
                </FormGroup>
              </FormControl>
            </div> : ''
          }

          <div className={classes.search}>
            <SearchBar classes={{autosuggestRoot: classes.searchBar}}
              fullWidth fullWidthInput={false} helperText={helperText}
              label="search"
              placeholder="enter atoms, codes, functionals, or other quantity values"
              data={data} searchValues={searchValues}
              InputLabelProps={{
                shrink: true
              }}
              onChanged={values => this.updateSearch({searchValues: values})}
            />
            <Divider className={classes.searchDivider} />
            <Tooltip title={showDetails ? 'hide statistics' : 'show statistics'}>
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

          <div className={classes.searchResults}>
            <Typography variant="caption" style={{margin: 12}}>
              About {total} results:
            </Typography>

            <SearchResultList
              data={data} total={total}
              onChange={this.updateSearchResultList}
              {...searchResultListState}
            />
          </div>
        </DisableOnLoading>
      </div>
    )
  }
}

export default compose(withApi(false), withErrors, withDomain, withStyles(SearchPage.styles))(SearchPage)
