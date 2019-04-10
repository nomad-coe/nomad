import React from 'react'
import PropTypes from 'prop-types'
import DFTSearchAggregations from './dft/DFTSearchAggregations'
import DFTEntryOverview from './dft/DFTEntryOverview'
import DFTEntryCards from './dft/DFTEntryCards'

const DomainContext = React.createContext()

export class DomainProvider extends React.Component {
  static propTypes = {
    children: PropTypes.oneOfType([
      PropTypes.arrayOf(PropTypes.node),
      PropTypes.node
    ]).isRequired
  }

  dft = {
    name: 'dft',
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
      total_energies: {
        label: 'Total energy calculations',
        renderResultString: count => (<span> with <b>{count}</b> total energy calculations</span>)
      },
      geometries: {
        label: 'Unique geometries',
        renderResultString: count => (<span> that simulate <b>{count}</b> unique geometries</span>)
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
  }

  state = {
    domain: this.dft
  }

  render() {
    return (
      <DomainContext.Provider value={this.state}>
        {this.props.children}
      </DomainContext.Provider>
    )
  }
}

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
