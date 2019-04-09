import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import Markdown from './Markdown'
import { kibanaBase, apiBase } from '../config'
import { compose } from 'recompose'
import { withApi } from './api'

class About extends React.Component {
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
          ## The nomad**@FAIR** prototype
          This is a prototype, a concept, for a continuation of the
          [NOMAD-coe](http://nomad-coe.eu) project. It is an attempt to redesign
          the nomad software and infrastructure with the following goals in mind:

          * more immediate and near-time use-modes (*staging area*, *code integration*, *on-site data*)
          * 3rd parties run instances of the nomad on their servers (*mirrors*, *oasis*, *industry* usage)
          * mirrors(partially) synchronize data with the central nomad instance (*data federation*)
          * the nomad architecture/infrastructure is used for related *domains* (e.g. experimental material science) or even more unrelated domains
          * nomad is integrated with existing *Open Data* initiatives and databases (*FAIRDI*, *EUDAT*, *optimade*)
          * we benefit from nomad being *Open Source* (public git, outside participation)

          ### Developer Documentation
          You find in depth developer documentation [here](${apiBase}/docs/index.html).
          It contains a general introduction to nomad, the underlying architecture,
          is (meta)data, and processing. Learn how to use the nomad ReST API. It
          contains information about how to develop nomad, how to operate it, how to
          contribute parser, and much more.

          ### ReST API
          Nomad services can also be accessed programmatically via nomad's
          ReST API. The API is described via [swagger](https://swagger.io/), therefore
          you can use your favorite swagger client library (e.g.
          [bravado](https://github.com/Yelp/bravado) for Python).
          Here is [our API's swagger UI](${apiBase}/) as reference documentation.

          ### Source code
          The source-code for nomad@FAIRDI is maintained at the MPCDF's
          [gitlab](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR). To push code, you
          need an MPCDF account and you can apply [here](https://www.mpcdf.mpg.de/userspace/forms/onlineregistrationform).

          ### Log management with Elastic stack
          We use a central logging system based on the *elastic*-stack
          (previously called *Elastic Logstash Kibana* (ELK)-stack).
          This system pushes logs, events, monitoring data,
          and other application metrics to a central database where it
          can be analysed visually. Here is the [link to Kibana](${kibanaBase}/)

          ### Test user
          During development this GUI might not be connected to the actual nomad
          repository. Therefore, you cannot create a user or login with an existing
          user. You might use our test users \`sheldon.cooper@nomad-fairdi.tests.de\`
          or \`leonard.hofstadter@nomad-fairdi.tests.de\` both
          with password \`password\`.

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

export default compose(withApi(), withStyles(About.styles))(About)
