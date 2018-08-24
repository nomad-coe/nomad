import React from 'react';
import PropTypes from 'prop-types';
import { withStyles, Paper, LinearProgress } from '@material-ui/core';
import ReactJson from 'react-json-view'
import api from '../api';
import Markdown from './Markdown';


class RepoCalc extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired
  }
  static styles = theme => ({
    root: {},
    calcData: {
      padding: theme.spacing.unit
    }
  });

  constructor(props) {
    super(props)
    this.state = {
      data: null
    }
  }

  componentDidMount() {
    const { uploadHash, calcHash} = this.props.match.params
    api.repo(uploadHash, calcHash).then(data => {
      this.setState({data: data})
    })
  }

  render() {
    const { classes } = this.props
    const { data } = this.state

    return (
      <div className={classes.root}>
        <Markdown>{`
          ## The Repository – Raw Code Data
        `}</Markdown>
        <Paper className={classes.calcData}>
          {
            data ?
              <ReactJson src={this.state.data} enableClipboard={false} collapsed={4} /> :
              <LinearProgress variant="query" />
          }
        </Paper>
      </div>

    )
  }
}

export default withStyles(RepoCalc.styles)(RepoCalc);