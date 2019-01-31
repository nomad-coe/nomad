import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core'
import Markdown from './Markdown'

class Home extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired
  }
  static styles = theme => ({
    root: {}
  });

  render() {
    const {classes} = this.props
    return (
      <div className={classes.root}>
        <Markdown>{`
          # The nomad**@FAIR** prototype
          This is a prototype, a concept, for a continuation of the
          [NOMAD-coe](http://nomad-coe.eu) project. It is an attempt to redesign
          the nomad software and infrastructure with the following goals in mind:

          * more immediate and near-time use-modes (*staging area*, *code integration*, *on-site data*)
          * 3rd parties run instances of the nomad on their servers (*mirrors*, *oasis*, *industry* usage)
          * mirrors(partially) synchronize data with the central nomad instance (*data federation*)
          * the nomad architecture/infrastructure is used for related *domains* (e.g. experimental material science) or even more unrelated domains
          * nomad is integrated with existing *Open Data* initiatives and databases (*FAIR*, *EUDAT*, *optimade*)
          * we benefit from nomad being *Open Source* (public git, outside participation)
        `}</Markdown>
      </div>

    )
  }
}

export default withStyles(Home.styles)(Home)
