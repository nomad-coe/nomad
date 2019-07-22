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
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired
  }

  domains = {
    DFT: {
      name: 'DFT',
      about: `
        ## The nomad**@FAIRDI** (*beta*)

        ### About nomad@FAIRDI

        After the conclusion of the original [NOMAD-coe](http://nomad-coe.eu) project,
        the newly founded NGO *FAIR Data Infrastructures* (FAIRDI) provides an
        umbrella to continue operation and further development of the Nomad
        material science data sharing platform.

        The immediate goal is to to consolidate and stabilize the nomad infrastructure, and
        as a first step, we refined the Nomad upload and data processing. This GUI introduces
        the *staging area* that allows you to observe your uploads processing and inspect
        the uploaded data before you decide to either publish your data or delete your upload again.

        Currently this is designed as a complement to the original [Nomad Repository GUI](https://repository.nomad-coe.eu/NomadRepository-1.1).
        You upload, process, inspect, and publish your data here. Here you have some
        capabilities to search and explore uploaded data. But to add comments, co-authors, and references,
        create data-sets, and manage your account you still have to use the original [Nomad Repository GUI](https://repository.nomad-coe.eu/NomadRepository-1.1).

        This GUI allows you to (menu on the left):
        * About: read about this, access the documentation, and API.
        * Search: inspect for existing data, your's and others (currently the search will only changed new data uploaded through nomad@FAIRDI).
        * Upload: drop data, view the processing, and publish your uploads.
        * Metainfo: browse the *metainfo*, Nomad's (meta-)data schema for processed data.

        ### How to use the new upload and processing

        **!Please read this, before you explore this new part of Nomad!**

        Feel free to explore this *new* part of Nomad, but expect that not everything will
        be working perfectly. Travel through the menu on the left and just
        use it. Feel free to upload data, look for limitations and things you do not like.
        The goal should be to figure out what is wrong and missing.

        Keep in mind that there are limitations:
        * You can only login with users that already exist in the Nomad Repository. If you
        are new to Nomad, visit the [Nomad Repository GUI](https://repository.nomad-coe.eu/NomadRepository-1.1)
        or register for a user account [here](http://nomad-repository.eu:8080/NomadRepository-1.1/register/).
        * When you published your data here, it will still take a day to index. Therefore,
        your data will not appear in the Nomad Repository immediately.
        * The existing entries from the original Nomad do not appear in the search. We
        are currently migrating the data. You will be able to see all your data, old and new, soon.

        For feedback and any issues you find, feel free to [open an issue](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/issues) or write
        an email to [markus.scheidgen@physik.hu-berlin.de](mailto:markus.scheidgen@physik.hu-berlin.de).
      `,
      entryLabel: 'calculation',
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
          renderResultString: count => (<span><b>{count.toLocaleString()}</b> entries</span>)
        },
        unique_code_runs: {
          label: 'Unique entries',
          tooltip: 'Counts duplicates only once.',
          renderResultString: count => (<span> and <b>{count.toLocaleString()}</b> unique entries</span>)
        },
        total_energies: {
          label: 'Total energy calculations',
          tooltip: 'Aggregates the number of total energy calculations as each entry can contain many calculations.',
          renderResultString: count => (<span> with <b>{count.toLocaleString()}</b> total energy calculations</span>)
        },
        geometries: {
          label: 'Unique geometries',
          tooltip: 'Aggregates the number of unique simulated system geometries in all entries.',
          renderResultString: count => (<span> that simulate <b>{count.toLocaleString()}</b> unique geometries</span>)
        },
        datasets: {
          label: 'Datasets',
          tooltip: 'Shows statistics in terms of datasets that entries belong to.',
          renderResultString: count => (<span> curated in <b>{count.toLocaleString()}</b> datasets</span>)
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
      about: `
        ## A Prototype for Experimental Material Science Data Sharing

        The original goal of the NOMAD CoE project was to provide a data sharing and
        publication platform for computational material science data. With this prototype,
        we want to apply Nomad ideas and implementations to experimental material science
        data.

        As a first step, this site demonstrates Nomad's \`domain specific\` search interface
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
