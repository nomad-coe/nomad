import React from 'react'
import PropTypes, { instanceOf } from 'prop-types'
import Markdown from '../Markdown'
import { withStyles, Paper, IconButton, FormGroup, Checkbox, FormControlLabel, FormLabel, Tooltip } from '@material-ui/core'
import UploadIcon from '@material-ui/icons/CloudUpload'
import Dropzone from 'react-dropzone'
import Upload from './Upload'
import { compose } from 'recompose'
import DeleteIcon from '@material-ui/icons/Delete'
import CheckIcon from '@material-ui/icons/Check'
import ReloadIcon from '@material-ui/icons/Cached'
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

To upload your own data, please put all the relevant files of all the calculations
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

When you publish your upload, it will take a night before it will appear in the
[NOMAD Repository](https://repository.nomad-coe.eu/NomadRepository-1.1/).
We are working on expediting this process.

#### Processing errors

We distinguish between uploads that fail processing completely and uploads that contain
entries that could not be processed. The former might be caused by issues during the
upload, bad file formats, etc. The latter (for more common) case means that not all of the provided
code input/output files could not be parsed by our parsers for various reasons.
The processing logs of the failed entries might provide some insight.

We do not allow the publishing of uploads that fail processing completely. Frankly, in most
cases there won't be any data to publish anyways. In the case of failed processing of
some entries, the data can still be published. You will be able to share it and create
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
      overflow: 'hidden'
    },
    formGroup: {
      paddingLeft: 0
    },
    uploadsLabel: {
      flexGrow: 1,
      paddingLeft: 0,
      padding: theme.spacing.unit * 2
    },
    uploads: {
      marginTop: theme.spacing.unit * 4
    },
    pagination: {
      textAlign: 'center'
    }
  })

  defaultData = {
    results: [],
    pagination: {
      total: 0,
      per_page: 10,
      page: 1
    }
  }

  state = {
    uploadCommand: {
      upload_command: 'loading ...',
      upload_tar_command: 'loading ...',
      upload_progress_command: 'loading ...'
    },
    data: {...this.defaultData},
    uploading: []
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

  update(newPage) {
    const { data: { pagination: { page, per_page }}} = this.state
    this.props.api.getUploads('all', newPage || page, per_page)
      .then(uploads => {
        this.setState({
          data: uploads,
          uploading: this.state.uploading.filter(upload => upload.current_task === 'uploading')})
      })
      .catch(error => {
        this.setState({data: {...this.defaultData}})
        this.props.raiseError(error)
      })
  }

  handleDoesNotExist(removedUpload) {
    const { uploading } = this.state
    this.setState({uploading: uploading.filter(upload => upload !== removedUpload)})
    this.update()
  }

  onDrop(files, rejectedFiles) {
    const upload = file => {
      const upload = this.props.api.createUpload(file.name)
      this.setState({uploading: [upload, ...this.state.uploading]})
      upload.uploadFile(file).catch(this.props.raiseError)
    }

    files.forEach(upload)
    rejectedFiles
      .filter(file => file.name.match(/(\.zip)|(\.bz)|(\.tgz)|(\.gz)|(\.bz2)$/i))
      .forEach(upload)
  }

  renderUploads() {
    const { classes } = this.props
    const { data: { results, pagination: { total, per_page, page }}, uploading } = this.state

    const renderUpload = upload => <Upload
      key={upload.gui_upload_id} upload={upload}
      onDoesNotExist={() => this.handleDoesNotExist(upload)}
    />

    return (<div className={classes.uploads}>
      <FormGroup className={classes.formGroup} row>
        <FormLabel className={classes.uploadsLabel}>Your uploads: </FormLabel>
        <Tooltip title="Reload uploads, e.g. after using the curl upload" >
          <IconButton onClick={() => this.update()}><ReloadIcon /></IconButton>
        </Tooltip>
      </FormGroup>
      {uploading.map(renderUpload)}
      {results.map(renderUpload)}
      {(total > per_page)
        ? <Pagination classes={{root: classes.pagination}}
          limit={per_page}
          offset={(page - 1) * per_page}
          total={total}
          onClick={(_, offset) => this.update((offset / per_page) + 1)}
          previousPageLabel={'prev'}
          nextPageLabel={'next'}
        /> : ''}
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
              'application/zip',
              'application/gzip',
              'application/bz2',
              'application/x-gzip',
              'application/x-bz2',
              'application/x-gtar',
              'application/x-tgz',
              'application/tar+gzip',
              'application/x-tar',
              'application/tar+bz2',
              'application/x-zip-compressed',
              'application/x-compressed',
              'application/x-zip']}
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
          <CopyToClipboard text={uploadCommand.upload_command} onCopy={() => null}>
            <Tooltip title="Copy command to clipboard">
              <IconButton>
                <ClipboardIcon />
              </IconButton>
            </Tooltip>
            {/* <button>Copy to clipboard with button</button> */}
          </CopyToClipboard>
          <HelpDialog icon={<MoreIcon/>} maxWidth="md" title="Alternative shell commands" content={`
            As an experienced shell and *curl* user, you can modify the commands to
            your liking.

            The given command can be modified. To see progress on large files, use
            \`\`\`
              ${uploadCommand.upload_progress_command}
            \`\`\`
            To \`tar\` and upload multiple folders in one command, use
            \`\`\`
            ${uploadCommand.upload_tar_command}
            \`\`\`

            ### Form data vs. streaming
            NOMAD accepts stream data (\`-T <local_file>\`) (like in the
            examples above) or multi-part form data (\`-X PUT -f file=@<local_file>\`):
            \`\`\`
            ${uploadCommand.upload_command_form}
            \`\`\`
            We generally recommend to use streaming, because form data can produce very
            large HTTP request on large files. Form data has the advantage of carrying
            more information (e.g. the file name) to our servers (see below).

            #### Upload names
            With multi-part form data (\`-X PUT -f file=@<local_file>\`), your upload will
            be named after the file by default. With stream data (\`-T <local_file>\`)
            there will be no default name. To set a custom name, you can use the URL
            parameter \`name\`:
            \`\`\`
            ${uploadCommand.upload_command_with_name}
            \`\`\`
            Make sure to user proper [URL encoding](https://www.w3schools.com/tags/ref_urlencode.asp)
            and shell encoding, if your name contains spaces or other special characters.
          `}/>
        </div>

        {this.renderUploads()}
      </div>
    )
  }
}

export default compose(withApi(true, false, 'To upload data, you must have a Nomad Repository account and you must be logged in.'), withCookies, withStyles(Uploads.styles))(Uploads)
