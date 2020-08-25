import React from 'react'
import DFTEntryOverview from './dft/DFTEntryOverview'
import DFTEntryCards from './dft/DFTEntryCards'
import EMSEntryOverview from './ems/EMSEntryOverview'
import EMSEntryCards from './ems/EMSEntryCards'
import {
  DFTSystemVisualizations, DFTPropertyVisualizations, DFTMethodVisualizations
} from './dft/DFTVisualizations'
import EMSVisualizations from './ems/EMSVisualizations'

/* eslint-disable react/display-name */

export const domains = ({
  dft: {
    name: 'computational',
    label: 'Computational material science data',
    key: 'dft',
    about: 'This include data from many computational material science codes',
    entryLabel: 'entry',
    entryLabelPlural: 'entries',
    entryTitle: data => data.dft && data.dft.code_name ? data.dft.code_name + ' run' : 'Code run',
    searchPlaceholder: 'enter atoms, codes, functionals, or other quantity values',
    /**
     * A set of components and metadata that is used to present tabs of search visualizations
     * in addition to the globally available elements and users view.
     */
    searchVisualizations: {
      'system': {
        component: DFTSystemVisualizations,
        label: 'System',
        description: 'Shows histograms on system metadata'
      },
      'method': {
        component: DFTMethodVisualizations,
        label: 'Method',
        description: 'Shows histograms on method metadata'
      },
      'properties': {
        component: DFTPropertyVisualizations,
        label: 'Properties',
        description: 'Shows histograms on the availability of key properties'
      }
    },
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
      'dft.calculations': {
        label: 'Single configuration calculations',
        shortLabel: 'SCC',
        tooltip: 'Aggregates the number of single configuration calculations (e.g. total energy calculations) as each entry can contain many calculations.',
        renderResultString: count => (<span> with <b>{count.toLocaleString()}</b> single configuration calculation{count === 1 ? '' : 's'}</span>)
      },
      // The unique_geometries search aggregates unique geometries based on 10^8 hashes.
      // This takes to long in elastic search for a reasonable user experience.
      // Therefore, we only support geometries without uniqueness check
      'dft.unique_geometries': {
        label: 'Unique geometries',
        shortLabel: 'Geometries',
        tooltip: 'Aggregates the number of simulated system geometries in all entries.',
        renderResultString: count => (<span> that simulate <b>{count.toLocaleString()}</b> unique geometrie{count === 1 ? '' : 's'}</span>)
      },
      'encyclopedia.material.materials': {
        label: 'Materials',
        tooltip: 'Shows statistics in terms of materials.',
        renderResultString: count => (<span> of <b>{count.toLocaleString()}</b> material{count === 1 ? '' : 's'}</span>)
      },
      datasets: {
        label: 'Datasets',
        tooltip: 'Shows statistics in terms of datasets that entries belong to.',
        renderResultString: count => (<span> curated in <b>{count.toLocaleString()}</b> dataset{count === 1 ? '' : 's'}</span>)
      }
    },
    defaultSearchMetric: 'code_runs',
    additionalSearchKeys: {
      raw_id: {},
      external_id: {},
      upload_id: {},
      calc_id: {},
      paths: {},
      pid: {},
      mainfile: {},
      calc_hash: {},
      formula: {},
      'dft.optimade': {},
      'dft.quantities': {},
      'dft.spacegroup': {},
      'dft.spacegroup_symbol': {},
      'dft.labels': {},
      upload_name: {}
    },
    /**
     * An dict where each object represents a column. Possible keys are label, render.
     * Default render
     */
    searchResultColumns: {
      'formula': {
        label: 'Formula',
        supportsSort: true
      },
      'dft.code_name': {
        label: 'Code',
        supportsSort: true
      },
      'dft.basis_set': {
        label: 'Basis set',
        supportsSort: true
      },
      'dft.xc_functional': {
        label: 'XC functionals',
        supportsSort: true
      },
      'dft.system': {
        label: 'System',
        supportsSort: true
      },
      'dft.crystal_system': {
        label: 'Crystal system',
        supportsSort: true
      },
      'dft.spacegroup_symbol': {
        label: 'Spacegroup',
        supportsSort: true
      },
      'dft.spacegroup': {
        label: 'Spacegroup (number)',
        supportsSort: true
      }
    },
    defaultSearchResultColumns: ['formula', 'dft.code_name', 'dft.system', 'dft.crystal_system', 'dft.spacegroup_symbol'],
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
    EntryCards: DFTEntryCards,
    /**
     * A component to render additional domain specific cards in the
     * the entry view. Needs to work with props: data (the entry data from the API),
     * loading (a bool with api loading status).
     */
    searchTabs: ['entries', 'materials', 'datasets', 'groups', 'uploads']
  },
  ems: {
    name: 'experimental',
    key: 'ems',
    label: 'Material science experiment data',
    about: 'This is metadata from material science experiments',
    entryLabel: 'entry',
    entryLabelPlural: 'entries',
    entryTitle: () => 'Experiment',
    searchPlaceholder: 'enter atoms, experimental methods, or other quantity values',
    searchVisualizations: {
      'metadata': {
        component: EMSVisualizations,
        label: 'Metadata',
        description: 'Shows histograms on system metadata'
      }
    },
    /**
     * Metrics are used to show values for aggregations. Each metric has a key (used
     * for API calls), a label (used in the select form), and result string (to show
     * the overall amount in search results).
     */
    searchMetrics: {
      code_runs: {
        label: 'Experiments',
        tooltip: 'Statistics will show the number of entires; usually each entry represents a single experiment.',
        renderResultString: count => (<span><b>{count}</b> entries</span>)
      },
      datasets: {
        label: 'Datasets',
        tooltip: 'Shows statistics in terms of datasets that entries belong to.',
        renderResultString: count => (<span> curated in <b>{count}</b> datasets</span>)
      }
    },
    defaultSearchMetric: 'code_runs',
    /**
     * An dict where each object represents a column. Possible keys are label, render.
     * Default render
     */
    searchResultColumns: {
      'formula': {
        label: 'Formula'
      },
      'ems.method': {
        label: 'Method'
      },
      'ems.experiment_location': {
        label: 'Location'
      },
      'ems.experiment_time': {
        label: 'Date/Time',
        render: entry => (entry.ems && entry.ems.experiment_time !== 'unavailable') ? new Date(entry.ems.experiment_time * 1000).toLocaleString() : 'unavailable'
      }
    },
    defaultSearchResultColumns: ['formula', 'ems.method', 'ems.experiment_location', 'ems.experiment_time'],
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
    EntryCards: EMSEntryCards,
    /**
     * Names of the possible search tabs for this domain
     */
    searchTabs: ['entries', 'datasets', 'uploads']
  }
})
