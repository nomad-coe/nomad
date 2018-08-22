import React from 'react';
import PropTypes from 'prop-types';
import { withStyles } from '@material-ui/core';


class Upload extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    upload: PropTypes.object.isRequired
  }
  static styles = theme => ({
      root: {
        background: 'red'
      }
  });

  constructor(props) {
    super(props)
    this.state = {
      upload: props.upload
    }
  }

  updateUpload() {
    window.setTimeout(() => {
      this.state.upload.update()
        .then(upload => {
          console.debug(`Sucessfully updated upload ${upload.id}.`)
          this.setState({upload: upload})
          if (!(upload.processing && upload.processing.status == 'SUCCESS')) {
            this.updateUpload()
          }
        })
    }, 1000)
  }

  componentDidMount() {
    this.updateUpload()
  }

  render() {
    const { classes } = this.props;
    const { upload } = this.state;

    return (
      <div className={classes.root}>
        <pre>
          {JSON.stringify(upload, null, 2)}
        </pre>
      </div>
    )
  }
}

export default withStyles(Upload.styles)(Upload);