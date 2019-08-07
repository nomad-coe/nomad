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
import MoreIcon from '@material-ui/icons/MoreHoriz'
import ClipboardIcon from '@material-ui/icons/Assignment'
import ConfirmDialog from './ConfirmDialog'
import HelpDialog from '../Help'
import { withApi } from '../api'
import { withCookies, Cookies } from 'react-cookie'
import Pagination from 'material-ui-flat-pagination'
import { CopyToClipboard } from 'react-copy-to-clipboard'

const publishedUploadsPageSize = 10

export const help = `
NOMAD now provides a two step upload process. After you upload your files, you
check NOMAD's processing of your files before you publish your data. This gives you
more control about how NOMAD will present your data.

#### Prepare and upload files

To upload your own data, please put all relevant files of all the calculations
you want to upload into a single \`*.zip\` or \`*.tar.gz\` archive.
We encourage you to add all code input and
output files, as well as any other auxiliary files that you might have created.
You can put data from multiple calculations, using your preferred directory
structure, into your archives. Drop your archive file(s) below. You can also
click the dropbox to select the file from your hard drive.

Alternatively, you can upload files via the given shell command.
Replace \`<local_file>\` with your archive file. After executing the command,
return here and press the reload button below).

There is a limit of 32 GB per upload. Please upload multiple archives, if
you have more than 32 GB of data to upload.

#### The staging area

Uploaded data will not be public immediately. We call this *staging area*.

Below you will find all your unpublished and published uploads.
The unpublished uploads are in the *staging area*. You can see the
progress on the processing, you can review your uploads, and publish or delete them again.

Click on an upload to see more details about its contents. Click on processed calculations
to see their metadata, archive data, and a processing log. Select uploads to *delete*
or *publish* them.

If you press publish, a dialog will appear that allows you to set an
*embargo* or publish your data as *Open Access* right away. The *embargo* allows you to shared
it with selected users, create a DOI for your data, and later publish the data.
The *embargo* might last up to 36 month before it becomes public automatically.
During an *embargo* some meta-data will be available.

When you published your upload, it will take a night before it will appear in the
[NOMAD Repository](https://repository.nomad-coe.eu/NomadRepository-1.1/).
We are working on improving this process.

#### Processing errors

We distinguish between uploads that fail processing completely and uploads that contain
entries that could not be processed. The former might be caused by issues during the
upload, bad file formats, etc. The latter (for more common) case means that not all of the provided
code input/output files could not be parsed by our parsers for various reasons.
The processing logs of the failed entries might provide some insight.

We do not allow to publish uploads that failed processing completely. Frankly, in most
cases there won't be any data to publish anyways. In the case of failed processing of
some entires, the data can still be published. You will be able to share it and create
DOIs for it, etc. The only shortcomings will be missing metadata (labeled *not processed*
or *unavailable*) and missing archive data. We continuously improve our parsers and
the missing information might be made available in the future.

#### Co-Authors, References, Comments, Datasets

Currently, this web-page is only about uploading your calculations. To further edit
comments, references, co-authors, share with other authors, or curate datasets, use
the [NOMAD Repository](https://repository.nomad-coe.eu/NomadRepository-1.1/) on
your published data (as usual).
`

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
    commandContainer: {
      display: 'flex',
      flexDirection: 'row',
      alignItems: 'center'
    },
    commandMarkup: {
      flexGrow: 1,
      marginRight: theme.spacing.unit,
      overflow: 'scroll'
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
    uploadsContainer: {
      marginTop: theme.spacing.unit * 4
    },
    pagination: {
      textAlign: 'center'
    }
  })

  state = {
    unpublishedUploads: null,
    publishedUploads: null,
    publishedUploadsPage: 1,
    publishedUploadsTotal: 0,
    uploadCommand: {
      upload_command: 'loading ...',
      upload_tar_command: 'loading ...',
      upload_progress_command: 'loading ...'
    },
    selectedUnpublishedUploads: [],
    showPublishDialog: false
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

  update(publishedUploadsPage) {
    this.props.api.getUnpublishedUploads()
      .then(uploads => {
        // const filteredUploads = uploads.filter(upload => !upload.is_state)
        this.setState({unpublishedUploads: uploads.results, selectedUnpublishedUploads: []})
      })
      .catch(error => {
        this.setState({unpublishedUploads: [], selectedUnpublishedUploads: []})
        this.props.raiseError(error)
      })
    this.props.api.getPublishedUploads(publishedUploadsPage, publishedUploadsPageSize)
      .then(uploads => {
        this.setState({
          publishedUploads: uploads.results,
          publishedUploadsTotal: uploads.pagination.total,
          publishedUploadsPage: uploads.pagination.page})
      })
      .catch(error => {
        this.setState({publishedUploads: []})
        this.props.raiseError(error)
      })
  }

  onDeleteClicked() {
    Promise.all(this.state.selectedUnpublishedUploads.map(upload => this.props.api.deleteUpload(upload.upload_id)))
      .then(() => this.update())
      .catch(error => {
        this.props.raiseError(error)
        this.update()
      })
  }

  onPublishClicked() {
    this.setState({showPublishDialog: true})
  }

  onPublish(withEmbargo) {
    Promise.all(this.state.selectedUnpublishedUploads
      .map(upload => this.props.api.publishUpload(upload.upload_id, withEmbargo)))
      .then(() => {
        this.setState({showPublishDialog: false})
        return this.update()
      })
      .catch(error => {
        this.props.raiseError(error)
        this.update()
      })
  }

  sortedUnpublishedUploads(order) {
    order = order || -1
    return this.state.unpublishedUploads.concat()
      .sort((a, b) => (a.gui_upload_id === b.gui_upload_id)
        ? 0
        : ((a.gui_upload_id < b.gui_upload_id) ? -1 : 1) * order)
  }

  handleDoesNotExist(nonExistingUpload) {
    this.setState({
      unpublishedUploads: this.state.unpublishedUploads.filter(upload => upload !== nonExistingUpload)
    })
  }

  handlePublished(publishedUpload) {
    this.update()
  }

  onDrop(files) {
    files.forEach(file => {
      const upload = this.props.api.createUpload(file.name)
      this.setState({unpublishedUploads: [...this.state.unpublishedUploads, upload]})
      upload.uploadFile(file).catch(this.props.raiseError)
    })
  }

  onSelectionChanged(upload, checked) {
    if (checked) {
      this.setState({selectedUnpublishedUploads: [upload, ...this.state.selectedUnpublishedUploads]})
    } else {
      const selectedUnpublishedUploads = [...this.state.selectedUnpublishedUploads]
      selectedUnpublishedUploads.splice(selectedUnpublishedUploads.indexOf(upload), 1)
      this.setState({selectedUnpublishedUploads: selectedUnpublishedUploads})
    }
  }

  onSelectionAllChanged(checked) {
    if (checked) {
      this.setState({selectedUnpublishedUploads: [...this.state.unpublishedUploads.filter(upload => !upload.tasks_running)]})
    } else {
      this.setState({selectedUnpublishedUploads: []})
    }
  }

  renderPublishedUploads() {
    const { classes } = this.props
    const { publishedUploadsTotal, publishedUploadsPage, publishedUploads } = this.state

    if (!publishedUploads || publishedUploads.length === 0) {
      return ''
    }

    return (<div className={classes.uploadsContainer}>
      <FormLabel className={classes.uploadsLabel}>Your published uploads: </FormLabel>
      <div className={classes.uploads}>
        <div>
          {
            publishedUploads.map(upload => (
              <Upload key={upload.gui_upload_id} upload={upload}
                checked={false}
                onCheckboxChanged={checked => true}/>
            ))
          }
          {
            (publishedUploadsTotal > publishedUploadsPageSize)
              ? <Pagination classes={{root: classes.pagination}}
                limit={publishedUploadsPageSize}
                offset={(publishedUploadsPage - 1) * publishedUploadsPageSize}
                total={publishedUploadsTotal}
                onClick={(_, offset) => this.update((offset / publishedUploadsPageSize) + 1)}
                previousPageLabel={'prev'}
                nextPageLabel={'next'}
              /> : ''
          }
        </div>
      </div>
    </div>)
  }

  renderUnpublishedUploads() {
    const { classes } = this.props
    const { selectedUnpublishedUploads, showPublishDialog } = this.state
    const unpublishedUploads = this.state.unpublishedUploads || []

    const reloadButton = <Tooltip title="Reload uploads, e.g. after using the curl upload" >
      <IconButton onClick={() => this.update()}><ReloadIcon /></IconButton>
    </Tooltip>

    return (<div className={classes.uploadsContainer}>
      <div style={{width: '100%'}}>
        {(unpublishedUploads.length === 0) ? ''
          : <FormLabel className={classes.uploadsLabel}>Your unpublished uploads: </FormLabel>
        }
        <FormGroup className={classes.selectFormGroup} style={{alignItems: 'center'}}row>
          {(unpublishedUploads.length === 0) ? <FormLabel label="all" style={{flexGrow: 1}}>You have currently no unpublished uploads</FormLabel>
            : <FormControlLabel label="all" style={{flexGrow: 1}} control={(
              <Checkbox
                checked={selectedUnpublishedUploads.length === unpublishedUploads.length && unpublishedUploads.length !== 0}
                onChange={(_, checked) => this.onSelectionAllChanged(checked)}
              />
            )} />
          }
          {reloadButton}
          <FormLabel classes={{root: classes.selectLabel}}>
            {`selected uploads ${selectedUnpublishedUploads.length}/${unpublishedUploads.length}`}
          </FormLabel>
          <Tooltip title="Delete selected uploads" >
            <div>
              <IconButton
                disabled={selectedUnpublishedUploads.length === 0}
                onClick={this.onDeleteClicked.bind(this)}
              >
                <DeleteIcon />
              </IconButton>
            </div>
          </Tooltip>

          <Tooltip title="Publish selected uploads" >
            <div>
              <IconButton
                disabled={selectedUnpublishedUploads.length === 0 || selectedUnpublishedUploads.some(upload => upload.tasks_status !== 'SUCCESS' || upload.total_calcs === 0)}
                onClick={() => this.onPublishClicked()}>
                <CheckIcon />
              </IconButton>
            </div>
          </Tooltip>

          <ConfirmDialog
            open={showPublishDialog}
            onClose={() => this.setState({showPublishDialog: false})}
            onPublish={(withEmbargo) => this.onPublish(withEmbargo)}
          />

        </FormGroup>
      </div>
      {
        (unpublishedUploads.length === 0)
          ? ''
          : <div className={classes.uploads}>
            { this.sortedUnpublishedUploads().map(upload => (
              <Upload key={upload.gui_upload_id} upload={upload}
                checked={selectedUnpublishedUploads.indexOf(upload) !== -1}
                onDoesNotExist={() => this.handleDoesNotExist(upload)}
                onPublished={() => this.handlePublished(upload)}
                onCheckboxChanged={checked => this.onSelectionChanged(upload, checked)}/>
            ))
            }
          </div>
      }
    </div>)
  }

  render() {
    const { classes } = this.props
    const { uploadCommand } = this.state

    return (
      <div className={classes.root}>
        <Paper className={classes.dropzoneContainer}>
          <Dropzone
            accept={[
              'application/zip', 'application/gzip', 'application/bz2', 'application/x-gzip',
              'application/x-bz2', 'application/x-gtar', 'application/x-tgz', 'application/tar+gzip',
              'application/x-tar', 'application/tar+bz2']}
            className={classes.dropzone}
            activeClassName={classes.dropzoneAccept}
            rejectClassName={classes.dropzoneReject}
            onDrop={this.onDrop.bind(this)}
          >
            <p>drop .tar.gz or .zip files here</p>
            <UploadIcon style={{fontSize: 36}}/>
          </Dropzone>
        </Paper>

        <div className={classes.commandContainer}>
          <div className={classes.commandMarkup}>
            <Markdown>{`
              \`\`\`
                ${uploadCommand.upload_command}
              \`\`\`
            `}</Markdown>
          </div>
          <CopyToClipboard text={uploadCommand} onCopy={() => null}>
            <Tooltip title="Copy command to clipboard">
              <IconButton>
                <ClipboardIcon />
              </IconButton>
            </Tooltip>
            {/* <button>Copy to clipboard with button</button> */}
          </CopyToClipboard>
          <HelpDialog icon={<MoreIcon/>} maxWidth="md" title="Alternative shell commands" content={`
            The given command can be modified. To see progress on large files, use
            \`\`\`
              ${uploadCommand.upload_progress_command}
            \`\`\`
            To \`tar\` and upload multiple folders in one command, use
            \`\`\`
            ${uploadCommand.upload_tar_command}
            \`\`\`
            As an experienced shell and *curl* user, you can modify the commands to
            your liking. NOMAD accepts stream data (\`-H <local_file>\`) or multi-form data (\`-f file=@<local_file>\`).
          `}/>
        </div>

        {this.renderUnpublishedUploads()}
        {this.renderPublishedUploads()}
      </div>
    )
  }
}

export default compose(withApi(true, false, 'To upload data, you must have a Nomad Repository account and you must be logged in.'), withCookies, withStyles(Uploads.styles))(Uploads)
