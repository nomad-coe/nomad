import React from 'react'
import PropTypes from 'prop-types'
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
    info: PropTypes.object,
    raiseError: PropTypes.func.isRequired
  }

  domains = {
    DFT: {
      name: 'DFT',
      about: `
        ### The New Nomad Data Upload

        This web-page was created to complement to the original
        [NOMAD Repository GUI](https://repository.nomad-coe.eu/NomadRepository-1.1).
        You upload, process, inspect, and publish your data here. Here you will be able
        to search and explore uploaded data. However, to add comments, co-authors, and references,
        create data-sets, and manage your account, you will still have to use the original
        [NOMAD Repository GUI](https://repository.nomad-coe.eu/NomadRepository-1.1).

        In the future, this web-page will include more and more features of other NOMAD
        components as an effort to consolidate the various web applications from the
        NOMAD Repository, Archive, Metainfo, Encyclopedia, and Analytics Toolkit.

        ### Limitations
        You can only login with users that already exist in the NOMAD Repository. If you
        are new to NOMAD, visit the [NOMAD Repository GUI](https://repository.nomad-coe.eu/NomadRepository-1.1)
        or register for a user account [here](http://nomad-repository.eu:8080/NomadRepository-1.1/register/).

        When you have published your data with the new upload process (menu item on the
        left of this page) you will still have to wait about a day for indexing to
        take place and have your calculation show up in the Nomad Repository GUI.
        Therefore, you should not expect your data to appear in the NOMAD Repository immediately.

        We have migrated all data from the NOMAD Repository to this new system. However, not
        all of the data was successfully processed by the new and more powerful parsers.
        We will continue to improve the parsers to raise the quality of archive data overtime.
        For some entries, no archive data might be currently available and some metadata might
        still be missing when you are exploring Nomad data using the new search and data
        exploring capabilities (menu items on the left).
      `,
      entryLabel: 'code run',
      searchPlaceholder: 'enter atoms, codes, functionals, or other quantity values',
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
          tooltip: 'The statistics will show the number of database entry. Each set of input/output files that represents a code run is an entry.',
          renderResultString: count => (<span><b>{count.toLocaleString()}</b> entr{count === 1 ? 'y' : 'ies'}</span>)
        },
        unique_entries: {
          label: 'Unique entries',
          tooltip: 'Counts duplicates only once.',
          renderResultString: count => (<span> and <b>{count.toLocaleString()}</b> unique entr{count === 1 ? 'y' : 'ies'}</span>)
        },
        // total_energies: {
        //   label: 'Total energy calculations',
        //   tooltip: 'Aggregates the number of total energy calculations as each entry can contain many calculations.',
        //   renderResultString: count => (<span> with <b>{count.toLocaleString()}</b> total energy calculation{count === 1 ? '' : 's'}</span>)
        // },
        calculations: {
          label: 'Single configuration calculations',
          tooltip: 'Aggregates the number of single configuration calculations (e.g. total energy calculations) as each entry can contain many calculations.',
          renderResultString: count => (<span> with <b>{count.toLocaleString()}</b> single configuration calculation{count === 1 ? '' : 's'}</span>)
        },
        unique_geometries: {
          label: 'Unique geometries',
          tooltip: 'Aggregates the number of unique simulated system geometries in all entries.',
          renderResultString: count => (<span> that simulate <b>{count.toLocaleString()}</b> unique geometrie{count === 1 ? '' : 's'}</span>)
        },
        datasets: {
          label: 'Datasets',
          tooltip: 'Shows statistics in terms of datasets that entries belong to.',
          renderResultString: count => (<span> curated in <b>{count.toLocaleString()}</b> dataset{count === 1 ? '' : 's'}</span>)
        }
      },
      additionalSearchKeys: {
        raw_id: {},
        upload_id: {},
        calc_id: {},
        paths: {},
        external_id: {},
        pid: {},
        mainfile: {},
        calc_hash: {},
        formula: {}
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
      about: `
        ## A Prototype for Experimental Material Science Data Sharing

        The original goal of the NOMAD CoE project was to provide a data sharing and
        publication platform for computational material science data. With this prototype,
        we want to apply NOMAD ideas and implementations to experimental material science
        data.

        As a first step, this site demonstrates NOMAD's \`domain specific\` search interface
        and how experiment (meta-)data can be represented. We want to explore what
        meta-data exists for material experiments, what is necessary to provide meaningful
        search capabilities, how we can implement FAIR data sharing principles, and
        how can we establish a community process to integrate the various experimental
        methods and respective data.
      `,
      entryLabel: 'experiment',
      searchPlaceholder: 'enter atoms, experimental methods, or other quantity values',
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
          tooltip: 'Statistics will show the number of database entry. Usually each entry represents a single experiment.',
          renderResultString: count => (<span><b>{count}</b> entries</span>)
        },
        datasets: {
          label: 'Datasets',
          tooltip: 'Shows statistics in terms of datasets that entries belong to.',
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
          render: time => time !== 'unavailable' ? new Date(time * 1000).toLocaleString() : time
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

  render() {
    const { info } = this.props

    return (
      <DomainContext.Provider value={{domain: info ? this.domains[info.domain.name] : this.domains.DFT}}>
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
