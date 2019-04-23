import React from 'react'
import PropTypes, { func } from 'prop-types'
import DFTSearchAggregations from './dft/DFTSearchAggregations'
import DFTEntryOverview from './dft/DFTEntryOverview'
import DFTEntryCards from './dft/DFTEntryCards'
import EMSSearchAggregations from './ems/EMSSearchAggregations'
import EMSEntryOverview from './ems/EMSEntryOverview'
import EMSEntryCards from './ems/EMSEntryCards'
import { withApi } from './api'

const DomainContext = React.createContext()

class DomainProviderBase extends React.Component {
  static propTypes = {
    children: PropTypes.oneOfType([
      PropTypes.arrayOf(PropTypes.node),
      PropTypes.node
    ]).isRequired,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired
  }

  domains = {
    DFT: {
      name: 'DFT',
      /**
       * A component that is used to render the search aggregations. The components needs
       * to work with props: aggregations (the aggregation data from the api),
       * searchValues (currently selected search values), metric (the metric key to use),
       * onChange (callback to propagate searchValue changes).
       */
      SearchAggregations: DFTSearchAggregations,
      /**
       * Metrics are used to show values for aggregations. Each metric has a key (used
       * for API calls), a label (used in the select form), and result string (to show
       * the overall amount in search results).
       */
      searchMetrics: {
        code_runs: {
          label: 'Entries',
          renderResultString: count => (<span><b>{count}</b> entries</span>)
        },
        unique_code_runs: {
          label: 'Unique entries',
          renderResultString: count => (<span> and <b>{count}</b> unique entries</span>)
        },
        total_energies: {
          label: 'Total energy calculations',
          renderResultString: count => (<span> with <b>{count}</b> total energy calculations</span>)
        },
        geometries: {
          label: 'Unique geometries',
          renderResultString: count => (<span> that simulate <b>{count}</b> unique geometries</span>)
        },
        datasets: {
          label: 'Datasets',
          renderResultString: count => (<span> curated in <b>{count}</b> datasets</span>)
        }
      },
      /**
       * An dict where each object represents a column. Possible keys are label, render.
       * Default render
       */
      searchResultColumns: {
        formula: {
          label: 'Formula'
        },
        code_name: {
          label: 'Code'
        },
        basis_set: {
          label: 'Basis set'
        },
        xc_functional: {
          label: 'XT treatment'
        },
        system: {
          label: 'System'
        },
        crystal_system: {
          label: 'Crystal system'
        },
        spacegroup_symbol: {
          label: 'Spacegroup'
        }
      },
      /**
       * A component to render the domain specific quantities in the metadata card of
       * the entry view. Needs to work with props: data (the entry data from the API),
       * loading (a bool with api loading status).
       */
      EntryOverview: DFTEntryOverview,
      /**
       * A component to render additional domain specific cards in the
       * the entry view. Needs to work with props: data (the entry data from the API),
       * loading (a bool with api loading status).
       */
      EntryCards: DFTEntryCards
    },
    EMS: {
      name: 'EMS',
      /**
       * A component that is used to render the search aggregations. The components needs
       * to work with props: aggregations (the aggregation data from the api),
       * searchValues (currently selected search values), metric (the metric key to use),
       * onChange (callback to propagate searchValue changes).
       */
      SearchAggregations: EMSSearchAggregations,
      /**
       * Metrics are used to show values for aggregations. Each metric has a key (used
       * for API calls), a label (used in the select form), and result string (to show
       * the overall amount in search results).
       */
      searchMetrics: {
        code_runs: {
          label: 'Entries',
          renderResultString: count => (<span><b>{count}</b> entries</span>)
        },
        datasets: {
          label: 'Datasets',
          renderResultString: count => (<span> curated in <b>{count}</b> datasets</span>)
        }
      },
      /**
       * An dict where each object represents a column. Possible keys are label, render.
       * Default render
       */
      searchResultColumns: {
        formula: {
          label: 'Formula'
        },
        method: {
          label: 'Method'
        },
        experiment_location: {
          label: 'Location'
        },
        experiment_time: {
          label: 'Date/Time',
          render: time => new Date(time * 1000).toLocaleString()
        }
      },
      /**
       * A component to render the domain specific quantities in the metadata card of
       * the entry view. Needs to work with props: data (the entry data from the API),
       * loading (a bool with api loading status).
       */
      EntryOverview: EMSEntryOverview,
      /**
       * A component to render additional domain specific cards in the
       * the entry view. Needs to work with props: data (the entry data from the API),
       * loading (a bool with api loading status).
       */
      EntryCards: EMSEntryCards
    }
  }

  state = {
    domain: this.domains.DFT
  }

  componentDidMount() {
    this.props.api.getInfo().then(info => {
      this.setState({domain: this.domains[info.domain.name] || this.domains.DFT})
    }).catch(error => {
      this.props.raiseError(error)
    })
  }

  render() {
    return (
      <DomainContext.Provider value={this.state}>
        {this.props.children}
      </DomainContext.Provider>
    )
  }
}

export const DomainProvider = withApi(false, false)(DomainProviderBase)

export function withDomain(Component) {
  function DomainConsumer(props) {
    return (
      <DomainContext.Consumer>
        {state => <Component {...state} {...props}/>}
      </DomainContext.Consumer>
    )
  }

  return DomainConsumer
}
