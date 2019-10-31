import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Card, CardHeader, CardContent } from '@material-ui/core'
import RawFiles from '../entry/RawFiles'

class DFTEntryCards extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired
  }

  static styles = theme => ({
    root: {}
  })

  render() {
    const { classes, ...props } = this.props

    return (
      <Card className={classes.root}>
        <CardHeader title="Raw files" />
        <CardContent classes={{root: classes.cardContent}}>
          <RawFiles {...props} />
        </CardContent>
      </Card>
    )
  }
}

export default withStyles(DFTEntryCards.styles)(DFTEntryCards)
