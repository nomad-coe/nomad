import React from 'react'
import PropTypes from 'prop-types'
import FileSaver from 'file-saver'
import { withApi } from '../api'
import { compose } from 'recompose'
import { withErrors } from '../errors'
import { apiBase } from '../../config'
import { withStyles, Tooltip } from '@material-ui/core'

class Download extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    fileName: PropTypes.string,
    url: PropTypes.string.isRequired,
    component: PropTypes.any,
    children: PropTypes.oneOfType([
      PropTypes.arrayOf(PropTypes.node),
      PropTypes.node
    ]).isRequired,
    api: PropTypes.object.isRequired,
    user: PropTypes.object,
    disabled: PropTypes.bool,
    raiseError: PropTypes.func.isRequired,
    tooltip: PropTypes.string,
    color: PropTypes.string,
    size: PropTypes.string
  }

  static styles = theme => ({
    root: {}
  })

  state = {
    preparingDownload: false
  }

  async onDownloadClicked() {
    const {url, api, user, fileName, raiseError} = this.props
    let fullUrl = `${apiBase}/${url}`
    let downloadUrl = fullUrl
    if (user) {
      api.getSignatureToken()
        .catch(error => {
          this.setState({preparingDownload: false})
          raiseError(error)
        })
        .then(result => {
          if (fullUrl.startsWith('/')) {
            fullUrl = `${window.location.origin}${fullUrl}`
          }
          const downloadUrl = new URL(fullUrl)
          downloadUrl.searchParams.append('signature_token', result)
          FileSaver.saveAs(downloadUrl.href, fileName)
          this.setState({preparingDownload: false})
        })
    } else {
      FileSaver.saveAs(downloadUrl, fileName)
      this.setState({preparingDownload: false})
    }
  }

  render() {
    const {classes, component, children, disabled, color, size, tooltip} = this.props
    const {preparingDownload} = this.state

    const Component = component

    const button = (
      <Component className={classes.root}
        disabled={disabled || preparingDownload} color={color} size={size}
        onClick={() => this.onDownloadClicked()}
      >
        {children}
      </Component>
    )

    if (tooltip && !disabled && !preparingDownload) {
      return <Tooltip title={tooltip}>{button}</Tooltip>
    } else {
      return button
    }
  }
}

export default compose(withApi(false), withErrors, withStyles(Download.styles))(Download)
