import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import { compose } from 'recompose'
import { withApi } from './api'
import Search from './search/Search'
import SearchContext from './search/SearchContext'

export const help = `
This page allows you to **inspect** and **manage** you own data.
`

class UserdataPage extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired
  }

  static styles = theme => ({
  })

  render() {
    return (
      <div>
        <SearchContext ownerTypes={['user', 'staging']} initialQuery={{owner: 'user'}} >
          <Search resultTab="entries"/>
        </SearchContext>
      </div>
    )
  }
}

export default compose(withApi(true, false, 'To manage you data, you must log in.'), withStyles(UserdataPage.styles))(UserdataPage)
