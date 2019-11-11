import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Fab, Card, CardContent, CardActions, Button, Typography } from '@material-ui/core'
import { Link } from 'react-router-dom'
import ReactJson from 'react-json-view'
import { compose } from 'recompose'
import Markdown from '../Markdown'
import { withApi } from '../api'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import Download from './Download'
import { ValueAttributes, MetaAttribute } from '../metaInfoBrowser/ValueCard'
import ApiDialogButton from '../ApiDialogButton'

export const help = `
The NOMAD **archive** provides data and meta-data in a common hierarchical format based on
well-defined quantity definitions that we call *metainfo*. This representation
is independent from the raw data format and provides a homogenous data stock.

You can click the various quantity values to see the quantity definition. Similarly,
you can click section names to get more information. Browse the *metainfo* to
learn more about NOMAD's archive format [here](/metainfo).
`

class ArchiveEntryView extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    info: PropTypes.object,
    raiseError: PropTypes.func.isRequired,
    uploadId: PropTypes.string.isRequired,
    calcId: PropTypes.string.isRequired
  }

  static styles = theme => ({
    root: {},
    error: {
      marginTop: theme.spacing.unit * 2
    },
    metaInfo: {
      height: '20vh',
      overflowY: 'auto',
      marginTop: theme.spacing.unit * 3,
      marginBottom: theme.spacing.unit * 3,
      display: 'flex',
      flexDirection: 'column'
    },
    metaInfoContent: {
      flex: 1
    },
    metaInfoActions: {
      flexDirection: 'row-reverse'
    },
    metaInfoDescription: {
      margin: `${theme.spacing.unit}px 0`
    },
    data: {
      height: '65vh',
      overflowX: 'auto'
    },
    downloadFab: {
      zIndex: 1,
      right: 32,
      bottom: 32,
      position: 'fixed !important'
    }
  })

  static defaultState = {
    data: null,
    showMetaInfo: false,
    doesNotExist: false
  }

  state = {
    metaInfo: null,
    ...ArchiveEntryView.defaultState
  }

  constructor(props) {
    super(props)
    this.unmounted = false
  }

  componentWillUnmount() {
    this.unmounted = true
  }

  componentDidMount() {
    this.updateArchive()
    this.updateMetaInfo()
  }

  componentDidUpdate(prevProps) {
    if (prevProps.api !== this.props.api || prevProps.info !== this.props.info) {
      this.updateMetaInfo()
    }

    if (prevProps.api !== this.props.api ||
        prevProps.uploadId !== this.props.uploadId ||
        prevProps.calcId !== this.props.calcId) {
      this.setState({...ArchiveEntryView.defaultState})
      this.updateArchive()
    }
  }

  updateMetaInfo() {
    if (this.props.api && this.props.info && !this.state.metaInfo) {
      this.props.api.getMetaInfo(this.props.info.domain.metainfo.all_package).then(metaInfo => {
        if (!this.unmounted) {
          this.setState({metaInfo: metaInfo})
        }
      })
    }
  }

  updateArchive() {
    const {uploadId, calcId, api} = this.props
    api.archive(uploadId, calcId).then(data => {
      if (!this.unmounted) {
        this.setState({data: data})
      }
    }).catch(error => {
      if (!this.unmounted) {
        this.setState({data: null})
      }
      if (error.name === 'DoesNotExist') {
        this.setState({doesNotExist: true})
      } else {
        this.props.raiseError(error)
      }
    })
  }

  handleShowMetaInfo(selection, more) {
    if (selection.name === '_name') {
      this.setState({showMetaInfo: selection.value})
    } else {
      this.setState({showMetaInfo: selection.name})
    }
  }

  render() {
    const { classes, uploadId, calcId } = this.props
    const { data, showMetaInfo, metaInfo, doesNotExist } = this.state
    const metaInfoData = metaInfo ? metaInfo.get(showMetaInfo) : null

    if (doesNotExist) {
      return (
        <Typography className={classes.error}>
          No archive does exist for this entry. Either the archive was not generated due
          to parsing or other processing errors (check the log tab), or the entry it
          self does not exist.
        </Typography>
      )
    }

    return (
      <div className={classes.root}>
        <Card className={classes.metaInfo}>
          <CardContent className={classes.metaInfoContent}>{
            showMetaInfo && metaInfo
              ? metaInfoData
                ? <div>
                  <MetaAttribute label={'metainfo name'} value={metaInfoData.name} />
                  <Markdown classes={{root: classes.metaInfoDescription}}>{metaInfoData.description}</Markdown>
                  <ValueAttributes definition={metaInfoData} />
                </div>
                : <Markdown>This value has **no** *meta-info* attached to it.</Markdown>
              : <Markdown>Click a value to show its *meta-info*!</Markdown>
          }</CardContent>
          <CardActions className={classes.metaInfoActions}>
            { (showMetaInfo && metaInfo && metaInfoData)
              ? <Button color="primary" component={props => <Link to={`/metainfo/${metaInfoData.name}`} {...props} />}>
                Goto Metainfo Browser
              </Button>
              : <Button color="primary" disabled>
                Goto Metainfo Browser
              </Button>}
            <ApiDialogButton
              component={props => (
                <Button color="primary" disabled={!(showMetaInfo && metaInfo && metaInfoData)} {...props}>
                  Metainfo JSON
                </Button>
              )}
              data={metaInfoData ? metaInfoData.miJson : {}} title="Metainfo JSON"
            />
          </CardActions>
        </Card>
        <Card className={classes.data}>
          <CardContent>
            {
              data && typeof data !== 'string'
                ? <ReactJson
                  src={this.state.data}
                  enableClipboard={false}
                  collapsed={2}
                  displayObjectSize={false}
                  onSelect={this.handleShowMetaInfo.bind(this)} />
                : <div>{
                  data
                    ? <div>
                      <Typography>Archive data is not valid JSON. Displaying plain text instead.</Typography>
                      <pre>{data || ''}</pre>
                    </div>
                    : <Typography>loading ...</Typography>
                }</div>
            }
          </CardContent>
        </Card>

        <Download
          classes={{root: classes.downloadFab}} tooltip="download calculation archive"
          component={Fab} className={classes.downloadFab} color="primary" size="medium"
          url={`archive/${uploadId}/${calcId}`} fileName={`${calcId}.json`}
        >
          <DownloadIcon />
        </Download>
      </div>
    )
  }
}

export default compose(withApi(false, true), withStyles(ArchiveEntryView.styles))(ArchiveEntryView)
