import React from 'react'
import PropTypes from 'prop-types'
import Markdown from './Markdown'
import { withStyles, Paper, IconButton, FormGroup, Checkbox, FormControlLabel, FormLabel,
  LinearProgress } from '@material-ui/core'
import UploadIcon from '@material-ui/icons/CloudUpload'
import Dropzone from 'react-dropzone'
import api from '../api'
import Upload from './Upload'
import { withErrors } from './errors'
import { compose } from 'recompose'
import DeleteIcon from '@material-ui/icons/Delete'
import CheckIcon from '@material-ui/icons/Check'
import ConfirmDialog from './ConfirmDialog'

class Uploads extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {
      width: '100%'
    },
    dropzoneContainer: {
      height: 192
    },
    dropzone: {
      textAlign: 'center',
      color: theme.palette.grey[500],
      fontSize: 24,
      height: '100%',
      display: 'flex',
      flexDirection: 'column',
      justifyContent: 'center',
      '& p': {
        marginTop: 0,
        marginBottom: theme.spacing.unit * 1
      },
      '& svg': {
        marginLeft: 'auto',
        marginRight: 'auto'
      }
    },
    dropzoneAccept: {
      background: theme.palette.primary.main,
      color: theme.palette.common.white
    },
    dropzoneReject: {
      background: 'red !important',
      color: theme.palette.common.white
    },
    selectFormGroup: {
      paddingLeft: theme.spacing.unit * 3
    },
    selectLabel: {
      padding: theme.spacing.unit * 2
    },
    uploads: {
      marginTop: theme.spacing.unit * 2
    }
  })

  state = {
    uploads: null,
    uploadCommand: null,
    selectedUploads: [],
    loading: true,
    showAccept: false
  }

  componentDidMount() {
    this.update()
    api.getUploadCommand()
      .then(command => this.setState({uploadCommand: command}))
      .catch(error => {
        this.props.raiseError(error)
      })
  }

  update() {
    this.setState({loading: true})
    api.getUploads()
      .then(uploads => {
        const filteredUploads = uploads.filter(upload => !upload.is_state)
        this.setState({uploads: filteredUploads, selectedUploads: [], loading: false})
      })
      .catch(error => {
        this.setState({uploads: [], selectedUploads: [], loading: false})
        this.props.raiseError(error)
      })
  }

  onDeleteClicked() {
    this.setState({loading: true})
    Promise.all(this.state.selectedUploads.map(upload => api.deleteUpload(upload.upload_id)))
      .then(() => this.update())
      .catch(error => {
        this.props.raiseError(error)
        this.update()
      })
  }

  onAcceptClicked() {
    this.setState({showAccept: true})
  }

  handleAccept() {
    this.setState({loading: true})
    Promise.all(this.state.selectedUploads.map(upload => api.unstageUpload(upload.upload_id)))
      .then(() => {
        this.setState({showAccept: false})
        return this.update()
      })
      .catch(error => {
        this.props.raiseError(error)
        this.update()
      })
  }

  onDrop(files) {
    files.forEach(file => {
      const upload = api.createUpload(file.name)
      this.setState({uploads: [...this.state.uploads, upload]})
      upload.uploadFile(file).catch(this.props.raiseError)
    })
  }

  onSelectionChanged(upload, checked) {
    if (checked) {
      this.setState({selectedUploads: [upload, ...this.state.selectedUploads]})
    } else {
      const selectedUploads = [...this.state.selectedUploads]
      selectedUploads.splice(selectedUploads.indexOf(upload), 1)
      this.setState({selectedUploads: selectedUploads})
    }
  }

  onSelectionAllChanged(checked) {
    if (checked) {
      this.setState({selectedUploads: [...this.state.uploads.filter(upload => upload.completed)]})
    } else {
      this.setState({selectedUploads: []})
    }
  }

  renderUploads() {
    const { classes } = this.props
    const { uploads, selectedUploads } = this.state

    if (uploads && uploads.length > 0) {
      return (
        <div>
          <div style={{width: '100%'}}>
            <Markdown text={'These are the *existing* uploads:'} />
            <FormGroup className={classes.selectFormGroup} row>
              <FormControlLabel label="all" style={{flexGrow: 1}} control={(
                <Checkbox
                  checked={selectedUploads.length === uploads.length}
                  onChange={(_, checked) => this.onSelectionAllChanged(checked)}
                />
              )} />
              <FormLabel classes={{root: classes.selectLabel}}>
                {`selected uploads ${selectedUploads.length}/${uploads.length}`}
              </FormLabel>
              <IconButton
                disabled={selectedUploads.length === 0}
                onClick={this.onDeleteClicked.bind(this)}
              >
                <DeleteIcon />
              </IconButton>

              <IconButton disabled={selectedUploads.length === 0} onClick={this.onAcceptClicked.bind(this)}>
                <CheckIcon />
              </IconButton>
              <ConfirmDialog open={this.state.showAccept} onClose={() => this.setState({showAccept: false})} onOk={this.handleAccept.bind(this)}>
                If you agree the selected uploads will move out of your private staging area into the public nomad.
              </ConfirmDialog>

            </FormGroup>
          </div>
          <div className={classes.uploads}>
            {this.state.uploads.map((upload) => (
              <Upload key={upload.upload_id} upload={upload}
                checked={selectedUploads.indexOf(upload) !== -1}
                onCheckboxChanged={checked => this.onSelectionChanged(upload, checked)}/>
            ))}
          </div>
        </div>
      )
    } else {
      return ''
    }
  }

  render() {
    const { classes } = this.props
    const { uploadCommand } = this.state

    return (
      <div className={classes.root}>
        <Markdown>{`
          ## Upload your own data
          You can upload your own data. Have your code output ready in a popular archive
          format (e.g. \`*.zip\` or \`*.tar.gz\`).  Your upload can
          comprise the output of multiple runs, even of different codes. Don't worry, nomad
          will find it, just drop it below:`}
        </Markdown>

        <Paper className={classes.dropzoneContainer}>
          <Dropzone
            accept="application/zip"
            className={classes.dropzone}
            activeClassName={classes.dropzoneAccept}
            rejectClassName={classes.dropzoneReject}
            onDrop={this.onDrop.bind(this)}
          >
            <p>drop files here</p>
            <UploadIcon style={{fontSize: 36}}/>
          </Dropzone>
        </Paper>

        <Markdown>{`
          Alternatively, you can upload files via the following shell command.
          Replace \`<local_file>\` with your file. After executing the command,
          return here and reload.

          \`\`\`
            ${uploadCommand}
          \`\`\`
          `}
        </Markdown>

        {this.renderUploads()}
        {this.state.loading ? <LinearProgress/> : ''}
      </div>
    )
  }
}

export default compose(withErrors, withStyles(Uploads.styles))(Uploads)
