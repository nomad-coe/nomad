import React from 'react'
import PropTypes from 'prop-types'
import Grid from '@material-ui/core/Grid'
import { Quantity } from './QuantityHistogram'
import { withStyles } from '@material-ui/core'
import SearchContext from './SearchContext'
import { compose } from 'recompose'
import { withApi } from '../api'

class UploadersList extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired
  }

  static styles = theme => ({
    root: {
      marginTop: theme.spacing(2)
    }
  })

  static contextType = SearchContext.type

  componentDidMount() {
    const {setStatisticsToRefresh} = this.context
    setStatisticsToRefresh('uploader')
  }

  render() {
    const {state: {usedMetric}} = this.context

    return (
      <Grid>
        <Quantity quantity="uploader" title="Uploaders" scale={1} metric={usedMetric} />
      </Grid>
    )
  }
}
export default compose(withApi(false, false), withStyles(UploadersList.styles))(UploadersList)
