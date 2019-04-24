import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Fab, Card, CardContent, CardActions, Button } from '@material-ui/core'
import { Link } from 'react-router-dom'
import ReactJson from 'react-json-view'
import { compose } from 'recompose'
import Markdown from '../Markdown'
import { withApi } from '../api'
import DownloadIcon from '@material-ui/icons/CloudDownload'
import Download from './Download'
import { ValueAttributes, MetaAttribute } from '../metaInfoBrowser/ValueCard'
import ApiDialogButton from '../ApiDialogButton'

class ArchiveEntryView extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    uploadId: PropTypes.string.isRequired,
    calcId: PropTypes.string.isRequired
  }

  static styles = theme => ({
    root: {},
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
      height: '60vh',
      overflowY: 'auto'
    },
    downloadFab: {
      zIndex: 1,
      right: 32,
      bottom: 32,
      position: 'fixed !important'
    }
  });

  constructor(props) {
    super(props)
    this.state = {
      data: null,
      metaInfo: null,
      showMetaInfo: false
    }
    this.unmounted = false
  }

  componentWillUnmount() {
    this.unmounted = true
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
    const {uploadId, calcId, api} = this.props
    api.archive(uploadId, calcId).then(data => {
      this.setState({data: data})
    }).catch(error => {
      this.setState({data: null})
      this.props.raiseError(error)
    })

    api.getMetaInfo().then(metaInfo => {
      if (!this.unmounted) {
        this.setState({metaInfo: metaInfo})
      }
    }).catch(error => {
      this.props.raiseError(error)
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
    const { data, showMetaInfo, metaInfo } = this.state
    const metaInfoData = metaInfo ? metaInfo.get(showMetaInfo) : null

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
              data
                ? <ReactJson
                  src={this.state.data}
                  enableClipboard={false}
                  collapsed={2}
                  displayObjectSize={false}
                  onSelect={this.handleShowMetaInfo.bind(this)} />
                : ''
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
