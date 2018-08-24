import React from 'react';
import PropTypes from 'prop-types';
import { withStyles, Paper, LinearProgress } from '@material-ui/core';
import ReactJson from 'react-json-view'
import api from '../api';
import Markdown from './Markdown';


class ArchiveCalc extends React.Component {
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
    api.archive(uploadHash, calcHash).then(data => {
      this.setState({data: data})
    })
  }

  render() {
    const { classes } = this.props
    const { data } = this.state

    return (
      <div className={classes.root}>
        <Markdown>{`
          ## The Archive – Code Independent Data
          The tree below shows all calculation data in nomad's *hierachical* and
          *code independent*. You can learn more about the different *sections* and
          *quantities* by visiting the [metainfo](/metainfo).
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

export default withStyles(ArchiveCalc.styles)(ArchiveCalc);