import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, Tab, Tabs } from '@material-ui/core'
import ArchiveEntryView from './ArchiveEntryView'
import ArchiveLogView from './ArchiveLogView'
import RepoEntryView from './RepoEntryView'
import { withApi, DoesNotExist } from '../api'
import { compose } from 'recompose'
import qs from 'qs'
import KeepState from '../KeepState'

class EntryPage extends React.Component {
  static styles = theme => ({
    root: {
    },
    content: {
      padding: `0 ${theme.spacing.unit * 3}px`,
      maxWidth: 1024,
      margin: 'auto'
    }
  })

  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    uploadId: PropTypes.string,
    calcId: PropTypes.string,
    location: PropTypes.object,
    query: PropTypes.bool
  }

  static defaultState = {
    viewIndex: 0,
    calcId: null,
    uploadId: null
  }

  state = {...EntryPage.defaultState}

  componentDidMount() {
    this.update()
  }

  componentDidUpdate(prevProps) {
    if (prevProps.query !== this.props.query ||
        prevProps.location.key !== this.props.location.key ||
        prevProps.uploadId !== this.props.uploadId ||
        prevProps.calcId !== this.props.calcId ||
        prevProps.query !== this.props.query ||
        prevProps.api !== this.props.api) {
      this.setState({ ...EntryPage.defaultState })
      this.update()
    }
  }

  update() {
    const { calcId, uploadId, query, location } = this.props
    if (query) {
      let queryParams = null
      if (location && location.search) {
        queryParams = qs.parse(location.search.substring(1))
      }
      this.props.api.search({...queryParams}).then(data => {
        if (data.results && data.results.length > 0) {
          const { calc_id, upload_id } = data.results[0]
          this.setState({uploadId: upload_id, calcId: calc_id, viewIndex: 0})
        } else {
          this.props.raiseError(new DoesNotExist())
        }
      }).catch(this.props.raiseError)
    } else {
      if (calcId && uploadId) {
        this.setState({calcId: calcId, uploadId: uploadId, viewIndex: 0})
      } else {
        // this should be unreachable
        this.props.raiseError(new DoesNotExist())
      }
    }
  }

  render() {
    const { classes } = this.props
    const { viewIndex, calcId, uploadId } = this.state

    if (calcId && uploadId) {
      const calcProps = { calcId: calcId, uploadId: uploadId }
      return (
        <div className={classes.root}>
          <Tabs
            className={classes.tabs}
            value={viewIndex}
            onChange={(event, state) => this.setState({viewIndex: state})}
            indicatorColor="primary"
            textColor="primary"
            variant="fullWidth"
          >
            <Tab label="Raw data" />
            <Tab label="Archive" />
            <Tab label="Logs" />
          </Tabs>

          <div className={classes.content}>
            <KeepState visible={viewIndex === 0} render={props => <RepoEntryView {...props} />} {...calcProps} />
            <KeepState visible={viewIndex === 1} render={props => <ArchiveEntryView {...props} />} {...calcProps} />
            <KeepState visible={viewIndex === 2} render={props => <ArchiveLogView {...props} />} {...calcProps} />
          </div>
        </div>
      )
    } else {
      return ''
    }
  }
}

export default compose(withApi(false, true), withStyles(EntryPage.styles))(EntryPage)
