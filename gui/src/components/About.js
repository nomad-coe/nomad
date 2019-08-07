import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import Markdown from './Markdown'
import { kibanaBase, apiBase, debug } from '../config'
import { compose } from 'recompose'
import { withApi } from './api'
import { withDomain } from './domains'

class About extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    domain: PropTypes.object.isRequired,
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
    const { classes, domain } = this.props
    const { info } = this.state

    return (
      <div className={classes.root}>
        <Markdown>{`
          ${domain.about}

          ### Getting Help
          If you encounter any difficulties, please write to
          [webmaster@nomad-repository.eu](mailto:webmaster@nomad-repository.eu). If you think
          that this web-page does not work as expected, or you want to start a discussion
          about possible features, feel free to open an issue on our [issue tracking
          system](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/issues).

          ### Developer Documentation
          You find in depth developer documentation [here](${apiBase}/docs/index.html).
          It contains a general introduction to NOMAD, the underlying architecture,
          is (meta)data, and processing. Learn how to use the NOMAD ReST API. It
          contains information about how to develop NOMAD, how to operate it, how to
          contribute parser, and much more.

          ### ReST API
          NOMAD services can also be accessed programmatically via NOMAD's
          ReST API. The API is described via [swagger](https://swagger.io/), therefore
          you can use your favorite swagger client library (e.g.
          [bravado](https://github.com/Yelp/bravado) for Python).
          Here is [our API's swagger UI](${apiBase}/) as reference documentation.

          ### Source code
          The source-code for this new version of NOMAD (dubbed *nomad@FAIRDI*) is maintained
          at the MPCDF's [gitlab](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR).
          To push code, you need an MPCDF account and you can apply
          [here](https://www.mpcdf.mpg.de/userspace/forms/onlineregistrationform).

          ${debug ? `
          ### Log management with Elastic stack
          We use a central logging system based on the *elastic*-stack
          (previously called *Elastic Logstash Kibana* (ELK)-stack).
          This system pushes logs, events, monitoring data,
          and other application metrics to a central database where it
          can be analysed visually. Here is the [link to Kibana](${kibanaBase}/)

          ### Test user
          During development this GUI might not be connected to the actual NOMAD
          repository. Therefore, you cannot create a user or login with an existing
          user. You might use the test user \`leonard.hofstadter@nomad-fairdi.tests.de\`
          with password \`password\`. The user \`sheldon.cooper@nomad-fairdi.tests.de\` is
          used for data that has no provenance with the original NOMAD CoE database.
          ` : ''}

          ### About this version
          - version: \`${info ? info.version : 'loading'}/${info ? info.release : 'loading'}\`
          - domain: ${info ? info.domain.name : 'loading'}
          - git: \`${info ? info.git.ref : 'loading'}; ${info ? info.git.version : 'loading'}\`
          - last commit message: *${info ? info.git.log : 'loading'}*
          - parsers: ${info ? info.parsers.join(', ') : 'loading'}
          - normalizers: ${info ? info.normalizers.join(', ') : 'loading'}
        `}</Markdown>
      </div>
    )
  }
}

export default compose(withApi(), withDomain, withStyles(About.styles))(About)
