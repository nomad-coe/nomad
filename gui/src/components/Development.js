import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import Markdown from './Markdown'
import { kibanaBase, apiBase } from '../config'
import { compose } from 'recompose'
import { withApi } from './api'

class Development extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {
      padding: theme.spacing.unit * 3
    }
  })

  state = {
    info: null
  }

  componentDidMount() {
    this.props.api.getInfo()
      .then(info => this.setState({info: info}))
      .catch(error => {
        this.props.raiseError(error)
      })
  }

  render() {
    const { classes } = this.props
    const { info } = this.state

    return (
      <div className={classes.root}>
        <Markdown>{`
          ### Nomad
          - version: \`${info ? info.version : 'loading'}/${info ? info.release : 'loading'}\`
          - domain: ${info ? info.domain.name : 'loading'}
          - git: \`${info ? info.git.ref : 'loading'}; ${info ? info.git.version : 'loading'}\`
          - last commit message: *${info ? info.git.log : 'loading'}*
          - parsers: ${info ? info.parsers.join(', ') : 'loading'}
          - normalizers: ${info ? info.normalizers.join(', ') : 'loading'}

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

export default compose(withApi(), withStyles(Development.styles))(Development)
