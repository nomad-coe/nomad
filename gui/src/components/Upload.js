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
      title: {
        minWidth: 200,
        marginRight: theme.spacing.unit * 2
      },
      stepper: {
        width: '100%',
        padding: 0
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
          console.debug(`Sucessfully updated upload ${upload.upload_id}.`)
          this.setState({upload: upload})
          if (upload.proc.status != 'SUCCESS') {
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

    const title = (
      <div className={classes.title}>
        <Typography variant="title">
          {upload.name || upload.upload_id}
        </Typography>
        <Typography variant="subheading">
          {new Date(Date.parse(upload.create_time)).toLocaleString()}
        </Typography>
      </div>
    );

    const batch = (
      <Typography className={classes.heading}>
        {upload.status}
      </Typography>
    )

    const proc = upload.proc
    console.assert(proc, 'Uploads always must have a proc')
    let activeStep = proc.task_names.indexOf(proc.current_task_name)
    if (proc.status == 'SUCCESS') {
      activeStep += 1
    }
    const stepper = (
      <Stepper activeStep={activeStep} classes={{root: classes.stepper}}>
        {proc.task_names.map((label, index) => (
          <Step key={label}>
            <StepLabel>{label}</StepLabel>
          </Step>
        ))}
      </Stepper>
    )

    return (
      <ExpansionPanel>
        <ExpansionPanelSummary expandIcon={<ExpandMoreIcon/>}>
          {title} {stepper}
        </ExpansionPanelSummary>
        <ExpansionPanelDetails style={{width: '100%'}}>
          <ReactJson src={upload} enableClipboard={false} collapsed={1}/>
        </ExpansionPanelDetails>
      </ExpansionPanel>
    )
  }
}

export default withStyles(Upload.styles)(Upload);