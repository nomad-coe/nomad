import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, LinearProgress } from '@material-ui/core'
import api from '../api'
import { compose } from 'recompose'
import { withErrors } from './errors'

class ArchiveLogView extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    uploadId: PropTypes.string.isRequired,
    calcId: PropTypes.string.isRequired
  }

  static styles = theme => ({
    root: {}
  });

  constructor(props) {
    super(props)
    this.state = {
      data: null
    }
  }

  componentDidMount() {
    const {uploadId, calcId} = this.props
    api.calcProcLog(uploadId, calcId).then(data => {
      this.setState({data: data})
    }).catch(error => {
      this.setState({data: null})
      this.props.raiseError(error)
    })
  }

  render() {
    const { classes } = this.props
    const { data } = this.state

    return (
      <div className={classes.root}>{
        data
          ? <pre>{data}</pre>
          : <LinearProgress variant="query" />
      }
      </div>
    )
  }
}

export default compose(withErrors, withStyles(ArchiveLogView.styles))(ArchiveLogView)
