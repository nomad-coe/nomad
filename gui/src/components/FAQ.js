import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import Markdown from './Markdown'

class FAQ extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired
  }

  static styles = theme => ({
    root: {
      padding: theme.spacing.unit * 3
    }
  })

  render() {
    const { classes } = this.props

    return (
      <div className={classes.root}>
        <Markdown>{`
          # Frequently Asked Questions (FAQ)

          ## Upload data, datasets, embargo, and DOIs

          ### How do I provide data and cite it in a paper?

          ### I don't want my data to be public yet, but I need to share it with other people.

          ### How can I share credit with my co authors?

          ### I uploaded data, but none of my files show up?

          ### Some of my data is marked as *not processed* or *unavailable*, what does this mean?

          ## Find and download data

          ## API use

          ## Contributing to NOMAD
        `}</Markdown>
      </div>
    )
  }
}

export default withStyles(FAQ.styles)(FAQ)
