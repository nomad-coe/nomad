import React from 'react'
import PropTypes from 'prop-types'
import { withApi } from '../api'
import { isEquivalent } from '../../utils'
import { compose } from 'recompose'
import { withDomain } from '../domains'

/**
 * A non visible component that keeps shared search state between all child components.
 */
class SearchContext extends React.Component {
  static propTypes = {
    query: PropTypes.object,
    initialQuery: PropTypes.object,
    ownerTypes: PropTypes.arrayOf(PropTypes.string)
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
    this.state.query = this.props.initialQuery || {}
  }

  defaultMetric = this.props.domain.defaultSearchMetric

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
    metric: this.defaultMetric,
    usedMetric: this.defaultMetric,
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
    const {api, raiseError, ownerTypes} = this.props
    const {request, query, metric} = this.state
    const search = {...request, ...query, metrics: metric === this.defaultMetric ? [] : [metric], ...(this.props.query || {})}

    api.search(search)
      .then(response => {
        // find the first statistic to determine which metric is used
        const {statistics} = response
        let usedMetric = this.defaultMetric
        const firstRealQuantitiy = Object.keys(statistics).find(key => key !== 'total')
        if (firstRealQuantitiy) {
          const firstValue = Object.keys(statistics[firstRealQuantitiy])[0]
          if (firstValue) {
            usedMetric = Object.keys(statistics[firstRealQuantitiy][firstValue])
              .find(metric => metric !== this.defaultMetric) || this.defaultMetric
          }
        }
        this.setState({response: response || SearchContext.emptyResponse, usedMetric: usedMetric})
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
        !isEquivalent(prevProps.query || {}, this.props.query || {})) {
      this.update()
    }
  }

  render() {
    const {children} = this.props
    const value = {
      state: this.state,
      props: this.props,
      setRequest: this.handleRequestChange,
      setQuery: this.handleQueryChange,
      setMetric: this.handleMetricChange
    }
    return <SearchContext.type.Provider value={value} >
      {children}
    </SearchContext.type.Provider>
  }
}

const withHoc = compose(withDomain, withApi(false, false))(SearchContext)
Object.assign(withHoc, {type: SearchContext.type})

export default withHoc
