
import React, { Component } from 'react'
import Viewer from './Viewer'
import PropTypes from 'prop-types'
import { withApi } from '../api'

class MetaInfoBrowser extends Component {
  static propTypes = {
    metainfo: PropTypes.string,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired
  }

  state = {
    metainfos: null
  }

  componentDidMount() {
    this.props.api.getMetaInfo().then(metainfos => {
      this.setState({metainfos: metainfos})
    }).catch(error => {
      this.props.raiseError(error)
    })
  }

  render() {
    const { metainfos } = this.state

    if (!metainfos) {
      return <div />
    }

    const metainfoName = this.props.metainfo || 'section_run'
    const metainfo = metainfos.resolve(metainfos.createProxy(metainfoName))

    return <Viewer rootElement={metainfo}
      packages={metainfos.contents}
      visiblePackages={metainfos.contents.map(pkg => pkg.name)} />
  }
}

export default withApi(false)(MetaInfoBrowser)
