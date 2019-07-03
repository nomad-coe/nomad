import React from 'react'
import PropTypes, { instanceOf } from 'prop-types'
import Markdown from '../Markdown'
import { withStyles, Paper, IconButton, FormGroup, Checkbox, FormControlLabel, FormLabel, Tooltip } from '@material-ui/core'
import UploadIcon from '@material-ui/icons/CloudUpload'
import Dropzone from 'react-dropzone'
import Upload from './Upload'
import { compose } from 'recompose'
import DeleteIcon from '@material-ui/icons/Delete'
import ReloadIcon from '@material-ui/icons/Cached'
import CheckIcon from '@material-ui/icons/Check'
import ConfirmDialog from './ConfirmDialog'
import { Help, Agree } from '../help'
import { withApi } from '../api'
import { withCookies, Cookies } from 'react-cookie'

class Uploads extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    cookies: instanceOf(Cookies).isRequired
  }

  static styles = theme => ({
    root: {
      padding: theme.spacing.unit * 3
    },
    dropzoneContainer: {
      height: 192,
      marginTop: theme.spacing.unit * 2,
      marginBottom: theme.spacing.unit * 2
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
    uploadCommand: 'loading ...',
    selectedUploads: [],
    showPublish: false
  }

  componentDidMount() {
    this.update()
    this.props.api.getUploadCommand()
      .then(command => {
        this.setState({uploadCommand: command})
      })
      .catch(error => {
        this.props.raiseError(error)
      })
  }

  update() {
    this.props.api.getUploads()
      .then(uploads => {
        const filteredUploads = uploads.filter(upload => !upload.is_state)
        this.setState({uploads: filteredUploads, selectedUploads: []})
      })
      .catch(error => {
        this.setState({uploads: [], selectedUploads: []})
        this.props.raiseError(error)
      })
  }

  onDeleteClicked() {
    Promise.all(this.state.selectedUploads.map(upload => this.props.api.deleteUpload(upload.upload_id)))
      .then(() => this.update())
      .catch(error => {
        this.props.raiseError(error)
        this.update()
      })
  }

  onPublishClicked() {
    this.setState({showPublish: true})
  }

  onPublish(withEmbargo) {
    Promise.all(this.state.selectedUploads
      .map(upload => this.props.api.publishUpload(upload.upload_id, withEmbargo)))
      .then(() => {
        this.setState({showPublish: false})
        return this.update()
      })
      .catch(error => {
        this.props.raiseError(error)
        this.update()
      })
  }

  sortedUploads(order) {
    order = order || -1
    return this.state.uploads.concat()
      .sort((a, b) => (a.gui_upload_id === b.gui_upload_id)
        ? 0
        : ((a.gui_upload_id < b.gui_upload_id) ? -1 : 1) * order)
  }

  handleDoesNotExist(nonExistingUupload) {
    this.setState({
      uploads: this.state.uploads.filter(upload => upload !== nonExistingUupload)
    })
  }

  onDrop(files) {
    files.forEach(file => {
      const upload = this.props.api.createUpload(file.name)
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
      this.setState({selectedUploads: [...this.state.uploads.filter(upload => !upload.tasks_running)]})
    } else {
      this.setState({selectedUploads: []})
    }
  }

  renderUploads() {
    const { classes } = this.props
    const { selectedUploads, showPublish } = this.state
    const uploads = this.state.uploads || []

    return (<div>
      <div style={{width: '100%'}}>
        <FormGroup className={classes.selectFormGroup} row>
          <FormControlLabel label="all" style={{flexGrow: 1}} control={(
            <Checkbox
              checked={selectedUploads.length === uploads.length && uploads.length !== 0}
              onChange={(_, checked) => this.onSelectionAllChanged(checked)}
            />
          )} />
          <Tooltip title="reload uploads" >
            <IconButton onClick={() => this.update()}><ReloadIcon /></IconButton>
          </Tooltip>
          <FormLabel classes={{root: classes.selectLabel}}>
            {`selected uploads ${selectedUploads.length}/${uploads.length}`}
          </FormLabel>
          <Tooltip title="delete selected uploads" >
            <div>
              <IconButton
                disabled={selectedUploads.length === 0}
                onClick={this.onDeleteClicked.bind(this)}
              >
                <DeleteIcon />
              </IconButton>
            </div>
          </Tooltip>

          <Tooltip title="publish selected uploads" >
            <div>
              <IconButton
                disabled={selectedUploads.length === 0 || selectedUploads.some(upload => upload.failed_calcs !== 0 || upload.total_calcs === 0)}
                onClick={() => this.onPublishClicked()}>
                <CheckIcon />
              </IconButton>
            </div>
          </Tooltip>

          <ConfirmDialog
            open={showPublish}
            onClose={() => this.setState({showPublish: false})}
            onPublish={(withEmbargo) => this.onPublish(withEmbargo)}
          />

        </FormGroup>
      </div>
      <div className={classes.uploads}>{
        (uploads.length > 0)
          ? (
            <div>
              <Help cookie="uploadList">{`
                These are all your uploads in the *staging area*. You can see the
                progress on data progresses and review your uploads before publishing
                them to the *nomad repository*.

                Select uploads to delete or publish them. Click on uploads to see individual
                calculations. Click on calculations to see more details on each calculation.

                When you select and click publish, you will be ask if you want to publish
                with or without the optional *embargo period*.
              `}</Help>
              {
                this.sortedUploads().map(upload => (
                  <Upload key={upload.gui_upload_id} upload={upload}
                    checked={selectedUploads.indexOf(upload) !== -1}
                    onDoesNotExist={() => this.handleDoesNotExist(upload)}
                    onCheckboxChanged={checked => this.onSelectionChanged(upload, checked)}/>
                ))
              }
            </div>
          ) : ''
      }</div>
    </div>)
  }

  render() {
    const { classes } = this.props
    const { uploadCommand } = this.state

    const agreement = `
      By uploading and downloading data, you agree to the
      [terms of use](https://www.nomad-coe.eu/the-project/nomad-repository/nomad-repository-terms).

      Note that uploaded files become downloadable by others. Uploaded data is licensed under the
      Creative Commons Attribution license ([CC BY 3.0](https://creativecommons.org/licenses/by/3.0/)).
      You can put an *embargo* on uploaded data. The *embargo period* lasts up to 36 month.
    `

    return (
      <div className={classes.root}>
        <Agree message={agreement} cookie="agreedToUploadTerms">
          <Help cookie="uploadHelp" component={Markdown}>{`
            To upload your own data, please put all relevant files of all the calculations
            you want to upload into a single \`*.zip\` or \`*.tar.gz\` archive.
            We encourage you to add all code input and
            output files, as well as any other auxiliary files that you might have created.
            You can put data from multiple calculations, using your preferred directory
            structure, into your archives. Drop your archive file(s) below. You can also
            click the dropbox to select the file from your hard drive.

            Uploaded data will not be public immediately. We call this *staging area*.
            After uploading and processing, you can decide if you want to make the data public,
            delete it again, or put an *embargo* on it.

            The *embargo* allows you to shared it with selected users, create a DOI
            for your data, and later publish the data. The *embargo* might last up to
            36 month before it becomes public automatically. During an *embargo*
            some meta-data will be available.

            There is a limit of 32 GB per upload. Please upload multiple archives, if
            you have more than 32 GB of data to upload.
          `}</Help>

          <Paper className={classes.dropzoneContainer}>
            <Dropzone
              accept={[
                'application/zip', 'application/gzip', 'application/bz2', 'application/x-gzip',
                'application/x-bz2', 'application/x-gtar', 'application/x-tgz', 'application/tar+gzip',
                'application/tar', 'application/tar+bz2']}
              className={classes.dropzone}
              activeClassName={classes.dropzoneAccept}
              rejectClassName={classes.dropzoneReject}
              onDrop={this.onDrop.bind(this)}
            >
              <p>drop .tar.gz or .zip files here</p>
              <UploadIcon style={{fontSize: 36}}/>
            </Dropzone>
          </Paper>

          <Help cookie="uploadCommandHelp">{`
            Alternatively, you can upload files via the following shell command.
            Replace \`<local_file>\` with your archive file. After executing the command,
            return here and press the reload button below). The same 32 GB limit applies.
          `}</Help>

          <Markdown>{`
            \`\`\`
              ${uploadCommand}
            \`\`\`
          `}</Markdown>

          {this.renderUploads()}
        </Agree>
      </div>
    )
  }
}

export default compose(withApi(true, false, 'To upload data, you must have a nomad account and you must be logged in.'), withCookies, withStyles(Uploads.styles))(Uploads)
