import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import { compose } from 'recompose'
import { withErrors } from '../errors'
import { withApi } from '../api'
import Search from './Search'
import SearchContext from './SearchContext'
import qs from 'qs'

export const help = `
This page allows you to **search** in NOMAD's data. The upper part of this page
gives you various options to enter and configure your search. The lower half
shows all data that fulfills your search criteria.

NOMAD's *domain-aware* search allows you to screen data by filtering based on
desired properties. This is different from basic *text-search* that traditional
search engines offer.

The search bar allows you to specify various quantity values that you want to
see in your results. This includes *authors*, *comments*, *atom labels*, *code name*,
*system type*, *crystal system*, *basis set types*, and *XC functionals*.

Alternatively, you can click *elements* or *metadata* to get a visual representation of
NOMAD's data as a periodic table or metadata charts. You can click the various
visualization elements to filter for respective quantities.

The visual representations show metrics for all data that fit your criteria.
You can display *entries* (e.g. code runs), *unique entries*, and *datasets*.
Other more specific metrics might be available.

The results tabs gives you a quick overview of all entries and datasets that fit your search.
You can click entries to see more details, download data, see the archive, etc. The *entries*
tab displays individual entries (e.g. code runs), the *grouped entries* tab will group
entries with similar metadata (e.g. it will group entries for the same material from the
  same user). The *dataset* tab, shows entry curated by user created datasets. You can
  click on datasets for a search page that will only display entries from the respective
  dataset.

The table columns can be configured. The *entries* tab also supports sorting. Selected
entries (or all entries) can be downloaded. The download will contain all user provided
raw calculation input and output files.

You can click entries to see more details about them. The details button will navigate
you to an entries page. These entry pages will show more metadata, raw files, the
entry's archive, and processing logs.
`

class SearchPage extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    user: PropTypes.object,
    location: PropTypes.object,
    raiseError: PropTypes.func.isRequired,
    update: PropTypes.number
  }

  static styles = theme => ({
    root: {
    },
    searchEntry: {
      padding: theme.spacing.unit * 3
    }
  })

  render() {
    const { classes, user, location, update } = this.props

    let query = {
      owner: 'all'
    }
    if (location && location.search) {
      query = {
        ...query,
        ...(qs.parse(location.search.substring(1)) || {})
      }
    }

    const withoutLogin = ['all']

    return (
      <div className={classes.root}>
        <SearchContext
          update={update}
          initialQuery={query}
          ownerTypes={['all', 'public'].filter(key => user || withoutLogin.indexOf(key) !== -1)}
        >
          <Search visualization="elements" tabs={['entries', 'groups', 'datasets']} />
        </SearchContext>
      </div>
    )
  }
}

export default compose(withApi(false), withErrors, withStyles(SearchPage.styles))(SearchPage)
