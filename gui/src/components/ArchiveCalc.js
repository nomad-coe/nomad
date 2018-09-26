import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Paper, LinearProgress, Typography } from '@material-ui/core'
import ReactJson from 'react-json-view'
import api from '../api'
import Markdown from './Markdown'
import { compose } from 'recompose'
import { withErrors } from './errors'

class ArchiveCalc extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    match: PropTypes.object.isRequired
  }
  static styles = theme => ({
    root: {},
    calcData: {
      padding: theme.spacing.unit
    },
    logs: {
      marginTop: theme.spacing.unit * 2,
      padding: theme.spacing.unit
    },
    metaInfo: {
      height: 120,
      padding: theme.spacing.unit * 2,
      overflowY: 'scroll'
    },
    metaInfoInstructions: {
      height: 100,
      paddingTop: 30,
      textAlign: 'center',
      color: 'grey'
    }
  });

  static metainfo = null

  constructor(props) {
    super(props)
    this.state = {
      data: null,
      logs: null,
      metaInfo: null,
      showMetaInfo: false
    }
  }

  componentDidMount() {
    const {uploadHash, calcHash} = this.props.match.params
    api.archive(uploadHash, calcHash).then(data => {
      this.setState({data: data})
    }).catch(error => {
      this.setState({data: null})
      this.props.raiseError(error)
    })

    api.calcProcLog(uploadHash, calcHash).then(logs => {
      if (logs && logs !== '') {
        this.setState({logs: logs})
      }
    })

    api.getMetaInfo().then(metaInfo => {
      this.setState({metaInfo: metaInfo})
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
    const { classes } = this.props
    const { data, showMetaInfo, metaInfo } = this.state
    const metaInfoData = metaInfo ? metaInfo[showMetaInfo] : null
    const { uploadHash, calcHash } = this.props.match.params
    return (
      <div className={classes.root}>
        <Markdown>{`
          ## The Archive – Code Independent Data
          All values in the archive data have a specific type and documentation
          associated with this type. This information is called the *meta-info*.
          You can learn more about the different *sections* and
          *quantities* by visiting the [meta-info](/metainfo) browser.
        `}</Markdown>
        <Paper className={classes.metaInfo}>
          {showMetaInfo && metaInfo
            ? metaInfoData
              ? <div>
                <Typography variant="title">{metaInfoData.name}</Typography>
                <Markdown>{metaInfoData.description}</Markdown>
              </div>
              : <div className={classes.metaInfoInstructions}>
                    this value has no meta-info attached to it
              </div>
            : <div className={classes.metaInfoInstructions}>
                click a value to show its meta-info
            </div>
          }
        </Paper>
        <Markdown>{`
          The tree below shows all calculation data in nomad's *hierachical* and
          *code independent* archive format. You can download it
          [here](${api.archiveUrl(uploadHash, calcHash)})
        `}</Markdown>
        <Paper className={classes.calcData}>
          {
            data
              ? <ReactJson
                src={this.state.data}
                enableClipboard={false}
                collapsed={4}
                displayObjectSize={false}
                onSelect={this.handleShowMetaInfo.bind(this)}/>
              : <LinearProgress variant="query" />
          }
        </Paper>
        { this.state.logs
          ?
          <Paper className={classes.logs}>
            <pre>
              {this.state.logs}
            </pre>
          </Paper>
          : ''}
      </div>

    )
  }
}

export default compose(withErrors, withStyles(ArchiveCalc.styles))(ArchiveCalc)
