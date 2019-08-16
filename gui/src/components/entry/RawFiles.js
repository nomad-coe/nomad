import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, FormGroup, FormControlLabel, Checkbox, FormLabel, IconButton, Divider, Typography } from '@material-ui/core'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import { withApi } from '../api'
import { compose } from 'recompose'
import Download from './Download'

class RawFiles extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    data: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    user: PropTypes.object,
    loading: PropTypes.number.isRequired,
    raiseError: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {},
    formLabel: {
      padding: theme.spacing.unit * 2
    }
  })

  state = {
    selectedFiles: [],
    uploadDirectory: null,
    files: null,
    doesNotExist: false
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
    const {data: {uploadId, calcId}} = this.props
    this.props.api.getRawFileListFromCalc(uploadId, calcId).then(data => {
      this.setState({files: data.contents, uploadDirectory: data.directory})
    }).catch(error => {
      this.setState({files: null})
      if (error.name === 'DoesNotExist') {
        this.setState({doesNotExist: true})
      } else {
        this.props.raiseError(error)
      }
    })
  }

  label(file) {
    return file
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
    const {classes, data: {upload_id, calc_id}, loading} = this.props
    const {selectedFiles, files, uploadDirectory, doesNotExist} = this.state

    const availableFiles = files ? files.map(file => file.name) : []

    const someSelected = selectedFiles.length > 0
    const allSelected = availableFiles.length === selectedFiles.length && someSelected

    if (doesNotExist) {
      return <Typography>
        The uploaded raw files for this entry do not exist. This is most likely a NOMAD
        issue. Please inform us, if this error persists.
      </Typography>
    }

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
            url={(selectedFiles.length === 1) ? `raw/${upload_id}/${uploadDirectory}/${selectedFiles[0]}` : `raw/${upload_id}?files=${encodeURIComponent(selectedFiles.map(file => `${uploadDirectory}/${file}`).join(','))}`}
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
