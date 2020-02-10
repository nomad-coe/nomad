import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, FormGroup, FormControlLabel, Checkbox, FormLabel, IconButton, Divider, Typography, Tooltip } from '@material-ui/core'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import { withApi } from '../api'
import { compose } from 'recompose'
import Download from './Download'
import ReloadIcon from '@material-ui/icons/Cached'
import ViewIcon from '@material-ui/icons/Search'
import InfiniteScroll from 'react-infinite-scroller'
import { ScrollContext } from '../App'

class RawFiles extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    uploadId: PropTypes.string.isRequired,
    calcId: PropTypes.string.isRequired,
    data: PropTypes.object,
    api: PropTypes.object.isRequired,
    user: PropTypes.object,
    loading: PropTypes.number.isRequired,
    raiseError: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {},
    formLabel: {
      padding: theme.spacing.unit * 2
    },
    shownFile: {
      color: theme.palette.secondary.main,
      overflowX: 'hidden',
      textOverflow: 'ellipsis',
      whiteSpace: 'nowrap',
      direction: 'rtl',
      textAlign: 'left'
    },
    fileContents: {
      width: '85%',
      overflowX: 'auto',
      color: theme.palette.primary.contrastText,
      backgroundColor: theme.palette.primary.dark,
      marginTop: theme.spacing.unit,
      padding: '3px 6px',
      fontFamily: 'Consolas, "Liberation Mono", Menlo, Courier, monospace',
      fontSize: 12
    },
    fileError: {
      marginTop: 16,
      padding: 8
    },
    fileNameFormGroup: {
      display: 'flex',
      flexWrap: 'nowrap'
    },
    fileNameFormGroupLabel: {
      flexGrow: 1,
      overflowX: 'hidden',
      marginRight: 0
    },
    fileNameLabel: {
      overflowX: 'hidden',
      textOverflow: 'ellipsis',
      whiteSpace: 'nowrap',
      direction: 'rtl',
      textAlign: 'left'
    }
  })

  static defaultState = {
    selectedFiles: [],
    fileContents: null,
    shownFile: null,
    files: null,
    doesNotExist: false
  }

  constructor(props) {
    super(props)
    this.handleFileClicked = this.handleFileClicked.bind(this)
  }

  state = {...RawFiles.defaultState}

  componentDidUpdate(prevProps) {
    if (prevProps.api !== this.props.api ||
        prevProps.uploadId !== this.props.uploadId ||
        prevProps.calcId !== this.props.calcId) {
      this.setState({...RawFiles.defaultState})
    }
  }

  update() {
    const { uploadId, calcId, raiseError } = this.props
    // this might accidentally happen, when the user logs out and the ids aren't
    // necessarily available anymore, but the component is still mounted
    if (!uploadId || !calcId) {
      return
    }

    this.props.api.getRawFileListFromCalc(uploadId, calcId).then(data => {
      const files = data.contents.map(file => data.directory ? `${data.directory}/${file.name}` : file.name)
      if (files.length > 500) {
        raiseError('There are more than 500 files in this entry. We can only show the first 500.')
      }
      this.setState({files: files})
    }).catch(error => {
      this.setState({files: null})
      if (error.name === 'DoesNotExist') {
        this.setState({doesNotExist: true})
      } else {
        raiseError(error)
      }
    })
  }

  label(file) {
    return file.split('/').reverse()[0]
  }

  handleSelectFile(file) {
    const {selectedFiles} = this.state
    const index = selectedFiles.indexOf(file)
    if (index === -1) {
      this.setState({selectedFiles: [file, ...selectedFiles]})
    } else {
      selectedFiles.splice(index, 1)
      this.setState({selectedFiles: selectedFiles})
    }
  }

  handleFileClicked(file) {
    const {api, uploadId, raiseError} = this.props
    this.setState({shownFile: file, fileContents: null})
    api.getRawFile(uploadId, file, {length: 16 * 1024})
      .then(contents => this.setState({fileContents: contents}))
      .catch(raiseError)
  }

  handleLoadMore(page) {
    const {api, uploadId, calcId, raiseError} = this.props
    const {fileContents, shownFile} = this.state

    // The infinite scroll component has the issue if calling load more whenever it
    // gets updates, therefore calling this infinitely before it gets any chances of
    // receiving the results (https://github.com/CassetteRocks/react-infinite-scroller/issues/163).
    // Therefore, we have to set hasMore to false first and set it to true again after
    // receiving actual results.
    this.setState({
      fileContents: {
        ...fileContents,
        hasMore: false
      }
    })

    if (fileContents.contents.length < (page + 1) * 16 * 1024) {
      api.getRawFile(uploadId, shownFile, {offset: page * 16 * 1024, length: 16 * 1024})
        .then(contents => {
          const {fileContents} = this.state
          // The back-button navigation might cause a scroll event, might cause to loadmore,
          // will set this state, after navigation back to this page, but potentially
          // different entry.
          if (this.props.calcId === calcId) {
            this.setState({
              fileContents: {
                ...contents,
                contents: ((fileContents && fileContents.contents) || '') + contents.contents
              }
            })
          }
        })
        .catch(error => {
          this.setState({fileContents: null, shownFile: null})
          raiseError(error)
        })
    }
  }

  filterPotcar(file) {
    if (file.toLowerCase().endsWith('potcar')) {
      return this.props.user && this.props.data.uploader.user_id === this.props.user.sub
    } else {
      return true
    }
  }

  render() {
    const {classes, uploadId, calcId, loading, data} = this.props
    const {selectedFiles, files, doesNotExist, fileContents, shownFile} = this.state

    const availableFiles = files || data.files || []

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
          {!files
            ? <Tooltip title="check for more files">
              <IconButton onClick={() => this.update()}>
                <ReloadIcon />
              </IconButton>
            </Tooltip> : ''
          }
          <FormLabel className={classes.formLabel}>
            {selectedFiles.length}/{availableFiles.length} files selected
          </FormLabel>
          <Download component={IconButton} disabled={selectedFiles.length === 0}
            color="secondary"
            tooltip="download selected files"
            url={(selectedFiles.length === 1) ? `raw/${uploadId}/${selectedFiles[0]}` : `raw/${uploadId}?files=${encodeURIComponent(selectedFiles.join(','))}&strip=true`}
            fileName={selectedFiles.length === 1 ? this.label(selectedFiles[0]) : `${calcId}.zip`}
          >
            <DownloadIcon />
          </Download>
        </FormGroup>
        <Divider />
        <div style={{display: 'flex', flexDirection: 'row'}}>
          <div style={{width: '25%'}}>
            {availableFiles.filter(this.filterPotcar.bind(this)).map((file, index) => (
              <FormGroup row key={index} className={classes.fileNameFormGroup}>
                <Tooltip title={file}>
                  <FormControlLabel
                    style={{flexGrow: 1, overflowX: 'hidden', textOverflow: 'ellipsis'}}
                    label={this.label(file)}
                    classes={{
                      root: classes.fileNameFormGroupLabel,
                      label: file === shownFile ? classes.shownFile : classes.fileNameLabel}}
                    control={
                      <Checkbox
                        disabled={loading > 0}
                        checked={selectedFiles.indexOf(file) !== -1}
                        onChange={() => this.handleSelectFile(file)} value={file}
                      />}
                  />
                </Tooltip>
                <Tooltip title='Show contents'>
                  <IconButton onClick={() => this.handleFileClicked(file)} color={file === shownFile ? 'secondary' : 'default'}>
                    <ViewIcon />
                  </IconButton>
                </Tooltip>
              </FormGroup>
            ))}
          </div>
          {fileContents && fileContents.contents !== null &&
            <ScrollContext.Consumer>
              {scroll =>
                <InfiniteScroll
                  className={classes.fileContents}
                  pageStart={0}
                  loadMore={this.handleLoadMore.bind(this)}
                  hasMore={fileContents.hasMore}
                  useWindow={false}
                  getScrollParent={() => scroll.scrollParentRef}
                >
                  <pre style={{margin: 0}}>
                    {`${fileContents.contents}`}
                    &nbsp;
                  </pre>
                </InfiniteScroll>
              }
            </ScrollContext.Consumer>
          }
          {fileContents && fileContents.contents === null &&
            <div className={classes.fileError}>
              <Typography color="error">
                Cannot display file due to unsupported file format.
              </Typography>
            </div>
          }
        </div>
      </div>
    )
  }
}

export default compose(withApi(false, true), withStyles(RawFiles.styles))(RawFiles)
