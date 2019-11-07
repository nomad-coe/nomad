import React from 'react'
import PropTypes from 'prop-types'
import { withApi } from '../api'
import { isEquivalent } from '../../utils'

/**
 * A non visible component that keeps shared search state between all child components.
 */
class SearchContext extends React.Component {
  static propTypes = {
    query: PropTypes.object
  }

  static emptyResponse = {
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

  static type = React.createContext()

  constructor(props) {
    super(props)
    this.handleRequestChange = this.handleRequestChange.bind(this)
    this.handleQueryChange = this.handleQueryChange.bind(this)
    this.handleMetricChange = this.handleMetricChange.bind(this)
  }

  state = {
    response: SearchContext.emptyResponse,
    request: {
      statistics: true,
      order_by: 'formula',
      order: 1,
      page: 1,
      per_page: 10,
      datasets: true,
      datasets_after: null
    },
    metric: 'code_runs',
    query: {}
  }

  handleRequestChange(changes) {
    this.setState({request: {...this.state.request, ...changes}})
  }

  handleQueryChange(changes, replace) {
    if (changes.atoms && changes.atoms.length === 0) {
      changes.atoms = undefined
    }
    if (changes.only_atoms && changes.only_atoms.length === 0) {
      changes.only_atoms = undefined
    }
    if (replace) {
      this.setState({query: changes})
    } else {
      this.setState({query: {...this.state.query, ...changes}})
    }
  }

  handleMetricChange(metric) {
    this.setState({metric: metric})
  }

  update() {
    const {api, raiseError} = this.props
    const {request, query, metric} = this.state
    const search = {...(this.props.query || {}), ...request, ...query, metrics: metric === 'code_runs' ? [] : [metric]}

    api.search(search)
      .then(response => {
        this.setState({response: response || SearchContext.emptyResponse})
      }).catch(error => {
        this.setState({response: SearchContext.emptyResponse})
        if (error.name !== 'NotAuthorized' || this.props.searchParameters.owner === 'all') {
          raiseError(error)
        }
      })
  }

  componentDidMount() {
    this.update()
  }

  componentDidUpdate(prevProps, prevState) {
    const {query, request, metric} = this.state
    if (
        prevState.query !== query ||
        prevState.request !== request ||
        prevState.metric !== metric ||
        !isEquivalent(prevProps.query, this.props.query)) {
      this.update()
    }
  }

  render() {
    const {children} = this.props
    const value = {
      state: this.state,
      setRequest: this.handleRequestChange,
      setQuery: this.handleQueryChange,
      setMetric: this.handleMetricChange
    }
    return <SearchContext.type.Provider value={value} >
      {children}
    </SearchContext.type.Provider>
  }
}

const withHoc = withApi(false, false)(SearchContext)
Object.assign(withHoc, {type: SearchContext.type})

export default withHoc
