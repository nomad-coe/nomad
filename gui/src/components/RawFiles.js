import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, FormGroup, FormControlLabel, Checkbox, FormLabel, IconButton, Divider } from '@material-ui/core'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import { withApi } from './api'
import { compose } from 'recompose'
import Download from './Download'

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

  render() {
    const {classes, files, uploadId, calcId} = this.props
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
          <Download component={IconButton} disabled={selectedFiles.length === 0}
            tooltip="download selected files"
            url={(selectedFiles.length === 1) ? `raw/${uploadId}/${selectedFiles[0]}` : `raw/${uploadId}?files=${encodeURIComponent(selectedFiles.join(','))}`}
            fileName={selectedFiles.length === 1 ? this.label(selectedFiles[0]) : `${calcId}.zip`}
          >
            <DownloadIcon />
          </Download>
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
