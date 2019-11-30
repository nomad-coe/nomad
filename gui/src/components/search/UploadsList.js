import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, TableCell, Toolbar, IconButton, FormGroup, Tooltip } from '@material-ui/core'
import { compose } from 'recompose'
import { withRouter } from 'react-router'
import { withDomain } from '../domains'
import NextIcon from '@material-ui/icons/ChevronRight'
import StartIcon from '@material-ui/icons/SkipPrevious'
import DataTable from '../DataTable'
import { withApi } from '../api'
import EditUserMetadataDialog from '../EditUserMetadataDialog'
import DownloadButton from '../DownloadButton'
import ClipboardIcon from '@material-ui/icons/Assignment'
import { CopyToClipboard } from 'react-copy-to-clipboard'
import DetailsIcon from '@material-ui/icons/MoreHoriz'

class UploadIdUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    uploadId: PropTypes.string.isRequired
  }

  static styles = theme => ({
    root: {
      display: 'inline-flex',
      alignItems: 'center',
      flexDirection: 'row',
      flexWrap: 'nowrap'
    }
  })

  render() {
    const {classes, uploadId} = this.props
    return <span className={classes.root}>
      {uploadId}
      <CopyToClipboard
        text={uploadId} onCopy={() => null}
      >
        <Tooltip title={`Copy to clipboard`}>
          <IconButton style={{margin: 3, marginRight: 0, padding: 4}}>
            <ClipboardIcon style={{fontSize: 16}} />
          </IconButton>
        </Tooltip>
      </CopyToClipboard>
    </span>
  }
}

export const UploadId = withStyles(UploadIdUnstyled.styles)(UploadIdUnstyled)

class UploadActionsUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    upload: PropTypes.object.isRequired,
    user: PropTypes.object,
    onChange: PropTypes.func,
    history: PropTypes.object.isRequired
  }

  static styles = theme => ({
    group: {
      flexWrap: 'nowrap',
      flexDirection: 'row-reverse'
    }
  })

  constructor(props) {
    super(props)
    this.handleEdit = this.handleEdit.bind(this)
    this.handleClickDetails = this.handleClickDetails.bind(this)
  }

  handleEdit() {
    const {onChange, upload} = this.props
    if (onChange) {
      onChange(upload)
    }
  }

  handleClickDetails() {
    this.props.history.push(`/uploads?open=${this.props.upload.example.upload_id}`)
  }

  render() {
    const {upload, user, classes} = this.props
    const editable = user && upload.example &&
      upload.example.authors.find(author => author.user_id === user.sub)

    const query = {upload_id: upload.example.upload_id}

    return <FormGroup row classes={{root: classes.group}}>
      <Tooltip title="Open this upload on the uploads page">
        <IconButton onClick={this.handleClickDetails}>
          <DetailsIcon />
        </IconButton>
      </Tooltip>
      {<DownloadButton query={query} tooltip="Download upload" />}
      {editable && <EditUserMetadataDialog
        title="Edit metadata of all entries in this upload"
        example={upload.example} query={query}
        total={upload.total} onEditComplete={this.handleEdit}
      />}
    </FormGroup>
  }
}

export const UploadActions = compose(withRouter, withApi(false), withStyles(UploadActionsUnstyled.styles))(UploadActionsUnstyled)

class UploadListUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    data: PropTypes.object,
    total: PropTypes.number,
    onChange: PropTypes.func.isRequired,
    history: PropTypes.any.isRequired,
    uploads_after: PropTypes.string,
    actions: PropTypes.element
  }

  static styles = theme => ({
    root: {
      overflow: 'auto',
      paddingLeft: theme.spacing.unit * 2,
      paddingRight: theme.spacing.unit * 2
    },
    scrollCell: {
      padding: 0
    },
    scrollBar: {
      minHeight: 56,
      padding: 0
    },
    scrollSpacer: {
      flexGrow: 1
    },
    clickableRow: {
      cursor: 'pointer'
    }
  })

  constructor(props) {
    super(props)
    this.renderEntryActions = this.renderEntryActions.bind(this)
  }

  columns = {
    upload_time: {
      label: 'Upload time',
      render: (upload) => new Date(upload.example.upload_time).toLocaleString()
    },
    upload_id: {
      label: 'Id',
      render: (upload) => <UploadId uploadId={upload.example.upload_id} />
    },
    last_processing: {
      label: 'Last processed',
      render: (upload) => new Date(upload.example.last_processing).toLocaleString()
    },
    version: {
      label: 'Processed with version',
      render: (upload) => upload.example.nomad_version
    },
    entries: {
      label: 'Entries',
      render: (upload) => upload.total
    },
    published: {
      label: 'Published',
      render: (upload) => upload.example.published ? 'Yes' : 'No'
    }
  }

  renderEntryActions(entry) {
    const {onChange} = this.props
    return <UploadActions search upload={entry} onChange={() => onChange({})} />
  }

  render() {
    const { classes, data, total, uploads_after, onChange, actions } = this.props
    const uploads = data.uploads || {values: []}
    const results = Object.keys(uploads.values).map(id => {
      return {
        id: id,
        total: uploads.values[id].total,
        example: uploads.values[id].examples[0]
      }
    })
    const per_page = 10
    const after = uploads.after

    let paginationText
    if (uploads_after) {
      paginationText = `next ${results.length} of ${total}`
    } else {
      paginationText = `1-${results.length} of ${total}`
    }

    const pagination = <TableCell colSpan={1000} classes={{root: classes.scrollCell}}>
      <Toolbar className={classes.scrollBar}>
        <span className={classes.scrollSpacer}>&nbsp;</span>
        <span>{paginationText}</span>
        <IconButton disabled={!uploads_after} onClick={() => onChange({uploads_after: null})}>
          <StartIcon />
        </IconButton>
        <IconButton disabled={results.length < per_page} onClick={() => onChange({uploads_after: after})}>
          <NextIcon />
        </IconButton>
      </Toolbar>
    </TableCell>

    return <DataTable
      entityLabels={['upload', 'uploads']}
      id={row => row.id}
      total={total}
      columns={this.columns}
      selectedColumns={['upload_time', 'upload_name', 'upload_id', 'entries', 'published']}
      entryActions={this.renderEntryActions}
      data={results}
      rows={per_page}
      actions={actions}
      pagination={pagination}
    />
  }
}

const UploadList = compose(withRouter, withDomain, withApi(false), withStyles(UploadListUnstyled.styles))(UploadListUnstyled)

export default UploadList
