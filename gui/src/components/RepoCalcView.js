import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, LinearProgress } from '@material-ui/core'
import ReactJson from 'react-json-view'
import api from '../api'
import { withErrors } from './errors'
import { compose } from 'recompose'

class RepoCalcView extends React.Component {
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
    api.repo(uploadId, calcId).then(data => {
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
          ? <ReactJson src={this.state.data} enableClipboard={false} collapsed={4} />
          : <LinearProgress variant="query" />
      }
      </div>
    )
  }
}

export default compose(withErrors, withStyles(RepoCalcView.styles))(RepoCalcView)
