import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import { FormControl, FormControlLabel, Checkbox, FormGroup, FormLabel, Tooltip } from '@material-ui/core'
import { compose } from 'recompose'
import { withErrors } from '../errors'
import { withApi, DisableOnLoading } from '../api'
import { appBase } from '../../config'
import Search from './Search';

export const help = `
This page allows you to **search** in NOMAD's data. The upper part of this page
gives you various options to enter and configure your search. The lower half
shows all data that fulfills your search criteria.

** Disclaimer: ** This is a preliminary version of the NOMAD software. It might
now show all of NOMAD's data. To see the full NOMAD dataset use the original
[NOMAD CoE Repository](https://repository.nomad-coe.eu/NomadRepository-1.1/search/)
for now.

#### Search Options

NOMAD's *domain-aware* search allows you to screen data by filtering based on
desired properties. This is different from basic *text-search* that traditional
search engines offer.

If you are logged-in, you can specify if you want to search among all data, publicly
available data, your own data, or just unpublished data in your [staging area](/uploads/).

The search bar allows you to specify various quantity values that you want to
see in your results. This includes *authors*, *comments*, *atom labels*, *code name*,
*system type*, *crystal system*, *basis set types*, and *XC functionals*.
Alternatively, you can click the periodic table and statistic bars to filter for respective
quantities.

The periodic table and bar-charts show metrics for all data that fit your criteria.
You can display *entries* (e.g. code runs), *unique entries*, and *datasets*.
Other more specific metrics might be available.

#### Search Results

The results table gives you a quick overview of all entries that fit your search.
You can click entries to see more details, download data, see the archive, etc.
The *raw files* tab, will show you all files that belong to the entry and offers a download
on individual, or all files. The *archive* tab, shows you the parsed data as a tree
data structure. This view is connected to NOMAD's [meta-info](${appBase}/metainfo), which acts a schema for
all parsed data. The *log* tab, will show you a log of the entry's processing.
`

class SearchPage extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    user: PropTypes.object,
    raiseError: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {
    },
    searchEntry: {
      padding: theme.spacing.unit * 3
    }
  })

  state = {
    owner: 'all'
  }

  render() {
    const { classes, user } = this.props
    const { owner } = this.state

    const ownerLabel = {
      all: 'All entries',
      public: 'Only public entries',
      user: 'Only your entries',
      staging: 'Staging area only'
    }

    const ownerTooltips = {
      all: 'This will show all entries in the database.',
      public: 'Do not show entries that are only visible to you.',
      user: 'Do only show entries visible to you.',
      staging: 'Will only show entries that you uploaded, but not yet published.'
    }

    const withoutLogin = ['all']

    return (
      <div className={classes.root}>
        <DisableOnLoading>
          <div className={classes.searchEntry}>
            <FormControl>
              <FormLabel>Filter entries and show: </FormLabel>
              <FormGroup row>
                {['all', 'public', 'user', 'staging']
                  .filter(key => user || withoutLogin.indexOf(key) !== -1)
                  .map(owner => (
                    <Tooltip key={owner} title={ownerTooltips[owner] + (user ? '' : 'You need to be logged-in for more options.')}>
                      <FormControlLabel
                        control={
                          <Checkbox checked={this.state.owner === owner} onChange={() => this.setState({owner: owner})} value="owner" />
                        }
                        label={ownerLabel[owner]}
                      />
                    </Tooltip>
                  ))}
              </FormGroup>
            </FormControl>
          </div>
        </DisableOnLoading>
        <Search searchParameters={{owner: owner}} showDetails />
      </div>
    )
  }
}

export default compose(withApi(false), withErrors, withStyles(SearchPage.styles))(SearchPage)
