import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import Markdown from './Markdown'
import gitInfo from '../gitinfo'
import { kibanaBase, apiBase } from '../config'

class Development extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired
  }

  static styles = theme => ({
    root: {}
  })

  render() {
    const { classes } = this.props

    return (
      <div className={classes.root}>
        <Markdown>{`
          ### Build info
          - version: \`${gitInfo.version}\`
          - ref: \`${gitInfo.ref}\`
          - last commit message: *${gitInfo.log}*

          ### ReST API
          Nomad services can also be accessed programatically via nomad's
          ReST API. The API is described via [swagger](https://swagger.io/), therefore
          you can use your favorit swagger client library (e.g.
          [bravado](https://github.com/Yelp/bravado) for Python).
          Here is [our API's swagger UI](${apiBase}/) as reference documentation.

          ### Elastic stack
          We use a central logging system based on the *elastic*-stack
          (previously called *Elastic Logstash Kibana* (ELK)-stack).
          This system pushes logs, events, monitoring data,
          and other application metrics to a central database where it
          can be analysed visually. Here is the [link to Kiaba](${kibanaBase}/)

          ### Test user
          During development this GUI might not be connected to the actual nomad
          repository. Therefore, you cannot create a user or login with an existing
          user. You might use our test users \`sheldon.cooper@nomad-fairdi.tests.de\`
          or \`leonard.hofstadter@nomad-fairdi.tests.de\` both
          with password \`password\`.
        `}</Markdown>
      </div>
    )
  }
}

export default withStyles(Development.styles)(Development)
