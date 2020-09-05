import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Card, CardHeader, CardContent } from '@material-ui/core'
import RawFiles from '../entry/RawFiles'
import Markdown from '../Markdown'

class EMSEntryCards extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    data: PropTypes.object.isRequired,
    loading: PropTypes.bool
  }

  static styles = theme => ({
    root: {},
    description: {
      marginBottom: theme.spacing(3)
    }
  })

  render() {
    const { classes, data, ...props } = this.props

    return (
      <Card className={classes.root}>
        <CardHeader title="Raw Data and Meta Data Files" />
        <CardContent classes={{root: classes.cardContent}}>
          <Markdown classes={{root: classes.description}}>{`
            The data for this experiment is externally stored and managed. Download the raw experiment data:
            [${data.ems && data.ems.repository_url}](${data.ems && data.ems.entry_repository_url}).

            The meta data describing this experiment in its original format, can be
            downloaded here directly:
          `}</Markdown>
          <RawFiles data={data} {...props} />
        </CardContent>
      </Card>
    )
  }
}

export default withStyles(EMSEntryCards.styles)(EMSEntryCards)
