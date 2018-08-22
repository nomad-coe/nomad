import React from 'react';
import PropTypes from 'prop-types';
import { withStyles, ExpansionPanel, ExpansionPanelSummary, Typography, ExpansionPanelDetails, Chip } from '@material-ui/core';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';


class Upload extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    upload: PropTypes.object.isRequired
  }
  static styles = theme => ({
      root: {},
      heading: {
        fontSize: theme.typography.pxToRem(15),
        fontWeight: theme.typography.fontWeightRegular,
      },
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
    }, 500)
  }

  componentDidMount() {
    this.updateUpload()
  }

  render() {
    const { classes } = this.props;
    const { upload } = this.state;

    const createTime = (
      <Typography className={classes.heading}>
        {new Date(Date.parse(upload.create_time)).toLocaleString()}
      </Typography>
    );

    const batch = (
      <Typography className={classes.heading}>
        {upload.status}
      </Typography>
    )

    return (
      <ExpansionPanel>
        <ExpansionPanelSummary expandIcon={<ExpandMoreIcon/>}>
          {createTime} {batch}
        </ExpansionPanelSummary>
        <ExpansionPanelDetails style={{width: '100%'}}>
          <Typography>
            <pre>
              {JSON.stringify(upload, null, 2)}
            </pre>
          </Typography>
        </ExpansionPanelDetails>
      </ExpansionPanel>
    )
  }
}

export default withStyles(Upload.styles)(Upload);