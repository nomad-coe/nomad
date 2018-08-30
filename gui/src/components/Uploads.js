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
import CommingSoon from './CommingSoon'

class Uploads extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {
      width: '100%'
    },
    dropzone: {
      textAlign: 'center',
      padding: theme.spacing.unit * 4,
      color: theme.palette.grey[500],
      fontSize: 24,
      '& p': {
        margin: theme.spacing.unit
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
    uploads: null, selectedUploads: [], loading: true, acceptCommingSoon: false
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
    this.setState({acceptCommingSoon: true})
  }

  onDrop(files) {
    files.forEach(file => {
      api.createUpload(file.name)
        .then(upload => upload.uploadFile(file))
        .then(upload => {
          this.setState({uploads: [...this.state.uploads, upload]})
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
      this.setState({selectedUploads: [...this.state.uploads]})
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
              <IconButton
                disabled={selectedUploads.length === 0}
                onClick={this.onAcceptClicked.bind(this)}
              >
                <CheckIcon />
              </IconButton>
              <CommingSoon
                open={this.state.acceptCommingSoon}
                onClose={() => this.setState({acceptCommingSoon: false})}>
                  This will allow you to accept uploads and their calculations.
                  Only accepted uploads will be available in nomad xt,
                  and they cannot be deleted anymore.
              </CommingSoon>
            </FormGroup>
          </div>
          <div className={classes.uploads}>
            {this.state.uploads.map((upload, key) => (
              <Upload key={key} upload={upload}
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

    return (
      <div className={classes.root}>
        <Markdown>{`
          ## Upload your own data to **nomad xt**

          You can upload your own data. Have your code output ready in a popular archive
          format (e.g. \`*.zip\` or \`*.tar.gz\`) and drop it below. Your upload can
          comprise the output of multiple runs, even of different codes. Don't worry, nomad
          will find it.
        `}</Markdown>
        <Paper>
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

        {this.renderUploads()}
        {this.state.loading ? <LinearProgress/> : ''}
      </div>
    )
  }
}

export default compose(withErrors, withStyles(Uploads.styles))(Uploads)
