import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import Markdown from './Markdown'
import { appBase, optimadeBase, apiBase, debug, consent } from '../config'
import { compose } from 'recompose'
import { withApi } from './api'
import packageJson from '../../package.json'
import { domains } from './domains'

class About extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    info: PropTypes.object,
    raiseError: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {
      padding: theme.spacing(3)
    }
  })

  render() {
    const { classes, info } = this.props

    return (
      <div className={classes.root}>
        <Markdown>{`
          # The NOMAD Repository and Archive

          This is the *graphical user interface* (GUI) for the NOMAD Repository and
          Archive. It allows you to **search, access, and download all NOMAD data** in its
          raw (Repository) and processed (Archive) form. You can **upload and manage your own
          raw materials science data**. Learn more about what data can be uploaded
          and how to prepare your data on the [NOMAD Repository homepage](https://repository.nomad-coe.eu/).
          You can access all published data without an account. If you want to provide
          your own data, please login or register for an account.

          In the future, this web-page will include more and more features of other NOMAD
          components as an effort to consolidate the various web applications from the
          NOMAD Repository, Archive, Metainfo, Encyclopedia, and Analytics Toolkit.


          ### Getting Help
          If you encounter any difficulties, please write to
          [webmaster@nomad-coe.eu](mailto:webmaster@nomad-coe.eu). If you think
          that this web-page is not working as expected, or if you want to start a discussion
          about possible features, feel free to open an issue on our [issue tracking
          system](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/issues).

          ### APIs
          The NOMAD Repository and Archive can also be accessed programmatically via ReST APIs.
          There is the proprietary NOMAD API and an implementation of the
          [OPTiMaDe API (0.10.0)](https://github.com/Materials-Consortia/OPTiMaDe/tree/master)
          standardized by the [OPTiMaDe consortium](https://www.optimade.org/)

          Both APIs are described via [swagger](https://swagger.io/) (also known as OpenAPI spec.),
          therefore you can use your favorite swagger client library
          (e.g. [bravado](https://github.com/Yelp/bravado) for Python).

          There are also web-based GUIs that allow to explore the APIs and their documentation:
          - [NOMAD API](${apiBase}/)
          - [OPTiMaDe API](${optimadeBase}/)

          There is a [tutorial on how to use the API with Python](${appBase}/docs/api_tutorial.html).

          There is also a Python library. You can use *pip* to install the library.
          \`\`\`
          pip install ${appBase}/dist/nomad-0.8.0.tar.gz
          \`\`\`

          The NOMAD Archive uses data that adheres to formal data definitions that we call
          the NOMAD Metainfo. You can download these definition in their [JSON form
          here](${apiBase}/archive/metainfo/all.nomadmetainfo.json). Otherwise, you
          can use the Meta Info browser from the menu to explore.

          ${debug ? `
          ### Material science data and domains
          Originally NOMAD was build for DFT calculations and data from the respective
          community code. By NOMAD supports multiple materials science domains:

          ${info && info.domains.map(domain => domains[domain.name]).map(domain => `- ${domain.name}: ${domain.about}`).join('\n')}
          ` : ''}

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

          ### Terms of use and licenses
          ${consent}

          ### About this version
          - version (API): \`${info ? info.version : 'loading'}/${info ? info.git.commit : 'loading'}\`
          - version (GUI): \`${packageJson.version}/${packageJson.commit}\`
          - domains: ${info ? Object.keys(info.domains).map(domain => info.domains[domain].name).join(', ') : 'loading'}
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

export default compose(withApi(), withStyles(About.styles))(About)
