import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, FormGroup, FormControlLabel, Checkbox, FormLabel, IconButton, Divider } from '@material-ui/core'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import { withApi } from '../api'
import { compose } from 'recompose'
import Download from './Download'

class RawFiles extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    uploadId: PropTypes.str.isRequired,
    mainfile: PropTypes.str.isRequired,
    api: PropTypes.object.isRequired,
    user: PropTypes.object,
    loading: PropTypes.number.isRequired
  }

  static styles = theme => ({
    root: {},
    formLabel: {
      padding: theme.spacing.unit * 2
    }
  })

  state = {
    selectedFiles: [],
    files: null
  }

  componentDidMount() {
    this.update()
  }

  componentDidUpdate(prevProps) {
    if (prevProps.api !== this.props.api) {
      this.update()
    }
  }

  update() {
    const {uploadId, mainfile} = this.props
    this.props.api.raw_file_list(uploadId, mainfile).then(data => {
      this.setState({files: data.files})
    }).catch(error => {
      this.setState({files: null})
      this.props.raiseError(error)
    })
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
    const {classes, uploadId, mainfile, loading} = this.props
    const {selectedFiles, files} = this.state

    const mainfileLocal = mainfile.split('/')[-1]

    let availableFiles = [
      {
        file: mainfileLocal,
        size: -1
      }
    ]

    if (files) {
      const mainfileIndex = files.findIndex(file => file.file === mainfileLocal)
    }

    const someSelected = selectedFiles.length > 0
    const allSelected = availableFiles.length === selectedFiles.length && someSelected

    return (
      <div className={classes.root}>
        <FormGroup row>
          <FormControlLabel
            label="select all" style={{flexGrow: 1}}
            control={
              <Checkbox value="select_all" checked={allSelected}
                indeterminate={!allSelected && someSelected}
                onChange={() => this.setState({selectedFiles: allSelected ? [] : availableFiles.slice()})}
              />
            }
          />
          <FormLabel className={classes.formLabel}>
            {selectedFiles.length}/{availableFiles.length} files selected
          </FormLabel>
          <Download component={IconButton} disabled={selectedFiles.length === 0}
            tooltip="download selected files"
            url={(selectedFiles.length === 1) ? `raw/${uploadId}/${selectedFiles[0]}` : `raw/${calc_id}?files=${encodeURIComponent(selectedFiles.join(','))}`}
            fileName={selectedFiles.length === 1 ? this.label(selectedFiles[0]) : `${calc_id}.zip`}
          >
            <DownloadIcon />
          </Download>
        </FormGroup>
        <Divider />
        <FormGroup row>
          {availableFiles.map((file, index) => (
            <FormControlLabel key={index} label={this.label(file)}
              control={
                <Checkbox
                  disabled={loading > 0}
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
