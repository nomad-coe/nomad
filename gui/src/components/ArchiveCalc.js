import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Paper, LinearProgress, Typography, Popover } from '@material-ui/core'
import ReactJson from 'react-json-view'
import api from '../api'
import Markdown from './Markdown'
import { compose } from 'recompose'
import { withErrors } from './errors'
import CalcProcLogPopper from './CalcProcLogPopper'

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
    },
    logLink: {
      fontSize: '1rem',
      lineHeight: '2',
      marginBlockStart: '1rem',
      marginBlockEnd: '1rem',
      '& a': {
        color: theme.palette.secondary.main,
        textDecoration: 'none',
        '&:hover': {
          textDecoration: 'underline'
        }
      }
    }
  });

  static metainfo = null

  state = {
    data: null,
    logs: null,
    metaInfo: null,
    showMetaInfo: false,
    showLogs: false
  }

  constructor(props) {
    super(props)
    this.logPopperAnchor = React.createRef()
  }

  componentDidMount() {
    const {uploadHash, calcHash} = this.props.match.params
    api.archive(uploadHash, calcHash).then(data => {
      this.setState({data: data})
    }).catch(error => {
      this.setState({data: null})
      this.props.raiseError(error)
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
      <div className={classes.root} ref={this.logPopperAnchor}>
        <Markdown>{`
          ## The Archive – Code Independent Data
          All values in the archive data have a specific type and documentation
          associated with this type. This information is called the *meta-info*.
          You can learn more about the different *sections* and
          *quantities* by visiting the [meta-info](/metainfo) browser.

          The tree below shows all calculation data in nomad's *hierachical* and
          *code independent* archive format. You can download it
          [here](${api.archiveUrl(uploadHash, calcHash)}). Click on values to
          see a *meta-info* description.
        `}</Markdown>
        <Typography className={classes.logLink}>
          The processing logs are available <a href="#" onClick={() => this.setState({showLogs: true})}>here</a>.
        </Typography>
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
        <CalcProcLogPopper
          open={this.state.showLogs}
          archiveId={`${uploadHash}/${calcHash}`}
          onClose={() => this.setState({showLogs: false})}
          anchorEl={this.logPopperAnchor.current}
          raiseError={this.props.raiseError}
        />
        <Popover
          open={(showMetaInfo && metaInfo && metaInfoData) ? true : false}
          anchorEl={this.logPopperAnchor.current}
          onClose={() => this.setState({showMetaInfo: null})}
          anchorOrigin={{
            vertical: 'center',
            horizontal: 'center',
          }}
          transformOrigin={{
            vertical: 'center',
            horizontal: 'center',
          }}
        >
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
        </Popover>
      </div>

    )
  }
}

export default compose(withErrors, withStyles(ArchiveCalc.styles))(ArchiveCalc)
