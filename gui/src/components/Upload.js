import React from 'react';
import PropTypes from 'prop-types';
import { withStyles, ExpansionPanel, ExpansionPanelSummary, Typography, ExpansionPanelDetails, Stepper, Step, StepLabel } from '@material-ui/core';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';
import ReactJson from 'react-json-view'


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

    const stepper = (
      <Stepper activeStep={upload.tasks.indexOf(upload.task)} orientation="vertical">
        {upload.tasks.map((label, index) => (
          <Step key={label}>
            <StepLabel>{label}</StepLabel>
          </Step>
        ))}
      </Stepper>
    )

    return (
      <ExpansionPanel>
        <ExpansionPanelSummary expandIcon={<ExpandMoreIcon/>}>
          {createTime} {stepper}
        </ExpansionPanelSummary>
        <ExpansionPanelDetails style={{width: '100%'}}>
          <ReactJson src={upload} enableClipboard={false} collapsed={1}/>
        </ExpansionPanelDetails>
      </ExpansionPanel>
    )
  }
}

export default withStyles(Upload.styles)(Upload);