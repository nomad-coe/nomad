import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, FormGroup, FormControlLabel, Checkbox, FormLabel, IconButton, Divider } from '@material-ui/core'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import FileSaver from 'file-saver'
import { apiBase } from '../config'
import { withApi } from './api'
import { compose } from 'recompose'

class RawFiles extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    uploadId: PropTypes.string.isRequired,
    calcId: PropTypes.string.isRequired,
    files: PropTypes.arrayOf(PropTypes.string).isRequired,
    api: PropTypes.object.isRequired,
    user: PropTypes.object
  }

  static styles = theme => ({
    root: {},
    formLabel: {
      padding: theme.spacing.unit * 2
    }
  })

  state = {
    selectedFiles: []
  }

  label(file) {
    return file.substring(file.lastIndexOf('/') + 1)
  }

  onSelectFile(file) {
    const {selectedFiles} = this.state
    const index = selectedFiles.indexOf(file)
    if (index === -1) {
      this.setState({selectedFiles: [file, ...selectedFiles]})
    } else {
      selectedFiles.splice(index, 1)
      this.setState({selectedFiles: selectedFiles})
    }
  }

  async onDownloadClicked() {
    const {uploadId, calcId, api, user} = this.props
    const files = this.state.selectedFiles
    const downloadFile = files.length === 1 ? this.label(files[0]) : `${calcId}.zip`

    let url
    let token
    if (user) {
      token = (await api.getSignatureToken()).token
      url = files.length === 1
        ? `${apiBase}/raw/${uploadId}/${files[0]}?token=${token}`
        : `${apiBase}/raw/${uploadId}?files=${encodeURIComponent(files.join(','))}&token=${token}`
    } else {
      url = files.length === 1
        ? `${apiBase}/raw/${uploadId}/${files[0]}`
        : `${apiBase}/raw/${uploadId}?files=${encodeURIComponent(files.join(','))}`
    }

    FileSaver.saveAs(url, downloadFile)
  }

  render() {
    const {classes, files} = this.props
    const {selectedFiles} = this.state
    const someSelected = selectedFiles.length > 0
    const allSelected = files.length === selectedFiles.length && someSelected

    return (
      <div className={classes.root}>
        <FormGroup row>
          <FormControlLabel
            label="select all" style={{flexGrow: 1}}
            control={
              <Checkbox value="select_all" checked={allSelected}
                indeterminate={!allSelected && someSelected}
                onChange={() => this.setState({selectedFiles: allSelected ? [] : files.slice()})}
              />
            }
          />
          <FormLabel className={classes.formLabel}>
            {selectedFiles.length}/{files.length} files selected
          </FormLabel>
          <IconButton
            disabled={selectedFiles.length === 0}
            onClick={() => this.onDownloadClicked()}
          >
            <DownloadIcon />
          </IconButton>
        </FormGroup>
        <Divider />
        <FormGroup row>
          {files.map((file, index) => (
            <FormControlLabel key={index} label={this.label(file)}
              control={
                <Checkbox
                  checked={selectedFiles.indexOf(file) !== -1}
                  onChange={() => this.onSelectFile(file)} value={file}
                />
              }
            />
          ))}
        </FormGroup>
      </div>
    )
  }
}

export default compose(withApi(false), withStyles(RawFiles.styles))(RawFiles)
