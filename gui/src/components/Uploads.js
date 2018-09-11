import React from 'react'
import PropTypes from 'prop-types'
import Markdown from './Markdown'
import { withStyles, Paper, IconButton, FormGroup, Checkbox, FormControlLabel, FormLabel,
  LinearProgress, InputLabel, Input, FormHelperText, Button, Popover, Grid, Typography,
  DialogContent, DialogActions} from '@material-ui/core'
import UploadIcon from '@material-ui/icons/CloudUpload'
import Dropzone from 'react-dropzone'
import api from '../api'
import Upload from './Upload'
import { withErrors } from './errors'
import { compose } from 'recompose'
import DeleteIcon from '@material-ui/icons/Delete'
import CheckIcon from '@material-ui/icons/Check'
import AddIcon from '@material-ui/icons/Add'
import UploadCommand from './UploadCommand'
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
    },
    uploadFormControl: {
      margin: theme.spacing.unit * 2
    },
    button: {
      margin: theme.spacing.unit
    },
    rightIcon: {
      marginLeft: theme.spacing.unit
    },
    uploadNameInput: {
      width: '100%'
    },
    uploadKindHeading: {
      paddingBottom: theme.spacing.unit
    },
    uploadKindDescription: {
      paddingTop: theme.spacing.unit,
      paddingBottom: theme.spacing.unit * 2
    },
    commandUpload: {
      height: 192
    }
  })

  state = {
    uploads: null,
    selectedUploads: [],
    loading: true,
    showAccept: false,
    uploadName: '',
    uploadCommand: null,
    showUploadCommand: false,
    uploadPopperAnchor: null
  }

  componentDidMount() {
    this.update()
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

  onCreateUploadCmdClicked(event) {
    event.persist()
    const existingUpload = this.state.uploads
      .find(upload => upload.name === this.state.uploadName && upload.waiting)
    if (existingUpload) {
      const upload = existingUpload
      this.setState({
        uploadCommand: upload.upload_command,
        showUploadCommand: true,
        uploadPopperAnchor: event.target})
    } else {
      api.createUpload(this.state.uploadName)
        .then(upload => {
          this.setState({
            uploads: [...this.state.uploads, upload],
            uploadCommand: upload.upload_command,
            showUploadCommand: true,
            uploadPopperAnchor: event.target})
        })
        .catch(error => {
          this.props.raiseError(error)
        })
    }
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
      api.createUpload(file.name)
        .then(upload => {
          this.setState({uploads: [...this.state.uploads, upload]})
          upload.uploadFile(file)
            .catch(this.props.raiseError)
        })
        .catch(this.props.raiseError)
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
    const { showUploadCommand, uploadCommand, uploadPopperAnchor } = this.state

    return (
      <div className={classes.root}>
        <Markdown>{`
          ## Upload your own data
          You can upload your own data. Have your code output ready in a popular archive
          format (e.g. \`*.zip\` or \`*.tar.gz\`).  Your upload can
          comprise the output of multiple runs, even of different codes. Don't worry, nomad
          will find it.`}
        </Markdown>
        <Grid container spacing={24}>
          <Grid item xs>
            <Typography variant="headline" className={classes.uploadKindHeading}>
              Browser upload
            </Typography>
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
            <Typography className={classes.uploadKindDescription}>
              Just drop your file above. You know this from many other services in the internet.
            </Typography>
          </Grid>
          <Grid item xs>
            <Typography variant="headline" className={classes.uploadKindHeading}>
              Command upload
            </Typography>
            <Paper className={classes.commandUpload}>
              <DialogContent>
                <InputLabel htmlFor="name-helper">Upload name</InputLabel>
                <Input className={classes.uploadNameInput}
                  id="name-helper" value={this.state.uploadName}
                  onChange={(event) => this.setState({uploadName: event.target.value})}
                />
                <FormHelperText id="name-helper-text">optiona, helps to track the upload</FormHelperText>
              </DialogContent>
              <DialogActions>
                <Button
                  color="primary" className={classes.button} variant="contained"
                  onClick={this.onCreateUploadCmdClicked.bind(this)}
                >
                    add upload
                  <AddIcon className={classes.rightIcon}/>
                </Button>
                <Popover
                  id="upload-command-popper"
                  onClose={() => this.setState({showUploadCommand: false})}
                  open={showUploadCommand}
                  anchorEl={uploadPopperAnchor}
                  anchorOrigin={{
                    vertical: 'bottom',
                    horizontal: 'center'
                  }}
                  transformOrigin={{
                    vertical: 'top',
                    horizontal: 'center'
                  }}
                >
                  <UploadCommand uploadCommand={uploadCommand} />
                </Popover>
              </DialogActions>
            </Paper>
            <Typography className={classes.uploadKindDescription}>
              You can upload your file via <strong>curl</strong>. Optionally, you
              can provide a name that will help to track different uploads.
              Without a name, you only have the upload time to follow your uploads.
              You can find the command by unfolding the new upload element.
            </Typography>
          </Grid>
        </Grid>

        {this.renderUploads()}
        {this.state.loading ? <LinearProgress/> : ''}
      </div>
    )
  }
}

export default compose(withErrors, withStyles(Uploads.styles))(Uploads)
