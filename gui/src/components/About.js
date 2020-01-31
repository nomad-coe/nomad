import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import Markdown from './Markdown'
import { appBase, optimadeBase, apiBase, debug, consent } from '../config'
import { compose } from 'recompose'
import { withApi } from './api'
import { withDomain } from './domains'
import packageJson from '../../package.json'

class About extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    info: PropTypes.object,
    domain: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {
      padding: theme.spacing.unit * 3
    }
  })

  render() {
    const { classes, domain, info } = this.props

    return (
      <div className={classes.root}>
        <Markdown>{`
          ${domain.about}

          ### Terms of use and licenses
          ${consent}

          ### Getting Help
          If you encounter any difficulties, please write to
          [webmaster@nomad-coe.eu](mailto:webmaster@nomad-coe.eu). If you think
          that this web-page is not working as expected, or if you want to start a discussion
          about possible features, feel free to open an issue on our [issue tracking
          system](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/issues).

          ### ReST APIs
          NOMAD services can also be accessed programmatically via ReST APIs.
          There is the proprietary NOMAD API and an implementation of the
          [OPTiMaDe API (0.10.0)](https://github.com/Materials-Consortia/OPTiMaDe/tree/master)
          standardized by the [OPTiMaDe consortium](https://www.optimade.org/)

          Both APIs are described via [swagger](https://swagger.io/) (also known as OpenAPI spec.),
          therefore you can use your favorite swagger client library
          (e.g. [bravado](https://github.com/Yelp/bravado) for Python).

          There are also web-based GUIs that allow to explore the APIs and their documentation:
          - [NOMAD API](${apiBase}/)
          - [OPTiMaDe API](${optimadeBase}/)

          ### Developer Documentation
          The [in-depth developer documentation](${appBase}/docs/index.html)
          contains a general introduction to NOMAD, the underlying architecture,
          is (meta)data, and processing. You will also find some information on how to use
          the NOMAD ReST API. It contains information about how to develop NOMAD, how to
          operate it, how to contribute parsers, and much more.

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
          can be analysed visually by us.

          ### Test user
          During development this GUI might not be connected to the actual NOMAD
          repository. Therefore, you cannot create a user or login with an existing
          user. You might use the test user \`leonard.hofstadter@nomad-fairdi.tests.de\`
          with password \`password\`. The user \`sheldon.cooper@nomad-fairdi.tests.de\` is
          used for data that has no provenance with the original NOMAD CoE database.
          ` : ''}

          ### About this version
          - version (API): \`${info ? info.version : 'loading'}/${info ? info.release : 'loading'}\`
          - version (GUI): \`${packageJson.version}\`
          - domain: ${info ? info.domain.name : 'loading'}
          - git: \`${info ? info.git.ref : 'loading'}; ${info ? info.git.version : 'loading'}\`
          - last commit message: *${info ? info.git.log : 'loading'}*
          - supported codes: ${info ? info.codes.join(', ') : 'loading'}
          - parsers: ${info ? info.parsers.join(', ') : 'loading'}
          - normalizers: ${info ? info.normalizers.join(', ') : 'loading'}
        `}</Markdown>
      </div>
    )
  }
}

export default compose(withApi(), withDomain, withStyles(About.styles))(About)
