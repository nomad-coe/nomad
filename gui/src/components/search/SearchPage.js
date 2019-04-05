import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import { FormControl, FormControlLabel, Checkbox, FormGroup,
  FormLabel, IconButton, MuiThemeProvider } from '@material-ui/core'
import { compose } from 'recompose'
import { withErrors } from '../errors'
import AnalyticsIcon from '@material-ui/icons/Settings'
import { analyticsTheme } from '../../config'
import Link from 'react-router-dom/Link'
import { withApi, DisableOnLoading } from '../api'
import SearchBar from './SearchBar'
import SearchResultList from './SearchResultList'
import SearchStatistics from './SearchStatistics'

class SearchPage extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    user: PropTypes.object,
    raiseError: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {
      padding: theme.spacing.unit * 3
    },
    selectFormGroup: {
      paddingLeft: theme.spacing.unit * 3
    },
    selectLabel: {
      padding: theme.spacing.unit * 2
    }
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
      ...SearchStatistics.defaultState
    },
    searchResultListState: {
      ...SearchResultList.defaultState
    }
  }

  constructor(props) {
    super(props)

    this.updateSearchResultList = this.updateSearchResultList.bind(this)
    this.updateSearch = this.updateSearch.bind(this)
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

  render() {
    const { classes, user } = this.props
    const { data, searchState, searchResultListState } = this.state
    const { searchValues } = searchState
    const { pagination: { total } } = data

    const ownerLabel = {
      all: 'All entries',
      public: 'Only public entries',
      user: 'Only your entries',
      staging: 'Only entries from your staging area'
    }

    return (
      <div className={classes.root}>
        <DisableOnLoading>
          { user
            ? <FormControl>
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
            </FormControl> : ''
          }

          <SearchBar
            fullWidth fullWidthInput={false} label="search" placeholder="enter atoms or other quantities"
            data={data} searchValues={searchValues}
            onChanged={values => this.updateSearch({searchValues: values})}
          />

          <SearchStatistics data={data} {...searchState} onChange={this.updateSearch} />

          <FormGroup className={classes.selectFormGroup} row>
            <FormLabel classes={{root: classes.selectLabel}} style={{flexGrow: 1}}>

            </FormLabel>
            <FormLabel classes={{root: classes.selectLabel}}>
            Analyse {total} code runs in an analytics notebook
            </FormLabel>
            <MuiThemeProvider theme={analyticsTheme}>
              <IconButton color="primary" component={Link} to={`/analytics`}>
                <AnalyticsIcon />
              </IconButton>
            </MuiThemeProvider>
          </FormGroup>

          <SearchResultList
            data={data} total={total}
            onChange={this.updateSearchResultList}
            {...searchResultListState}
          />
        </DisableOnLoading>
      </div>
    )
  }
}

export default compose(withApi(false), withErrors, withStyles(SearchPage.styles))(SearchPage)
