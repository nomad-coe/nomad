import React from 'react';
import Markdown from './Markdown';
import { withStyles, Paper } from '@material-ui/core';
import UploadIcon from '@material-ui/icons/CloudUpload';
import Dropzone from 'react-dropzone';
import api from '../api';
import Upload from './Upload'
import { withErrors } from './errors';
import { compose } from 'recompose';

var styles = theme => ({
    root: {
      width: '100%',
    },
    dropzone: {
      textAlign: 'center',
      padding: theme.spacing.unit * 4,
      color: theme.palette.grey[500],
      fontSize: 24,
      '& p': {
        margin: theme.spacing.unit
      }
    },
    dropzoneAccept: {
      background: theme.palette.primary.main,
      color: theme.palette.common.white
    },
    dropzoneReject: {
      background: 'red !important',
      color: theme.palette.common.white
    },
    uploads: {
      marginTop: theme.spacing.unit * 2
    }
});

class Uploads extends React.Component {

  constructor(props) {
    super(props)
    this.state = {
      uploads: []
    }
  }

  componentDidMount() {
    api.getUploads()
      .then(uploads => {
        const filteredUploads = uploads.filter(upload => !upload.is_state)
        this.setState({uploads: filteredUploads})
      })
      .catch(error => {
        this.setState({uploads: []})
        this.props.raiseError(error)
      })
  }

  onDrop(files) {
    files.forEach(file => {
      api.createUpload(file.name)
        .then(upload => upload.uploadFile(file))
        .then(upload => {
          this.setState({uploads: [...this.state.uploads, upload]})
        })
        .catch(this.props.raiseError)
    });
  }

  render() {
    const { classes } = this.props;

    return (
      <div className={classes.root}>
        <Markdown>{`
          ## Upload your own data to **nomad xt**

          You can upload your own data. Have your code output ready in a popular archive
          format (e.g. \`*.zip\` or \`*.tar.gz\`) and drop it below. Your upload can
          comprise the output of multiple runs, even of different codes. Don't worry, nomad
          will find it.
        `}</Markdown>
        <Paper>
          <Dropzone
              accept="application/zip"
              className={classes.dropzone}
              activeClassName={classes.dropzoneAccept}
              rejectClassName={classes.dropzoneReject}
              onDrop={this.onDrop.bind(this)}
          >
            <p>drop files here</p>
            <UploadIcon style={{fontSize: 36}}/>
          </Dropzone>
        </Paper>
        {this.state.uploads.length > 0 ?
          <Markdown text={'These are the *existing* uploads:'} />:
          ''
        }
        <div className={classes.uploads}>
          {this.state.uploads.map((upload, key) => (<Upload key={key} upload={upload}/>))}
        </div>
      </div>
    )
  }
}

export default compose(withErrors, withStyles(styles))(Uploads)