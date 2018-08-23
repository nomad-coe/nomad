import React from 'react';
import PropTypes from 'prop-types';
import { withStyles, ExpansionPanel, ExpansionPanelSummary, Typography, ExpansionPanelDetails, Stepper, Step, StepLabel, Table, TableRow, TableCell } from '@material-ui/core';
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
      details: {
        padding: 0,
        display: 'block',
      },
      detailsContent: {
        margin: theme.spacing.unit * 3
      },
      title: {
        flexBasis: '20%',
        flexShrink: 0,
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
        {proc.task_names.map((label, index) => {
          let optional = null;
          if (proc.task_names[index] === 'parse_all') {
            label = 'parse'
            if (proc.calc_procs.length > 0) {
              optional = (
                <Typography variant="caption">
                  {proc.calc_procs.filter(p => p.status === 'SUCCESS').length}/{proc.calc_procs.length}
                </Typography>
              );
            }
          }
          return (
            <Step key={label}>
              <StepLabel optional={optional}>{label}</StepLabel>
            </Step>
          )
        })}
      </Stepper>
    )

    return (
      <ExpansionPanel>
        <ExpansionPanelSummary expandIcon={<ExpandMoreIcon/>}>
          {title} {stepper}
        </ExpansionPanelSummary>
        <ExpansionPanelDetails style={{width: '100%'}} classes={{root: classes.details}}>
          {(proc.calc_procs.length === 0) ? <Typography className={classes.detailsContent}>No calculcations found.</Typography> : (
            <Table>
              {proc.calc_procs.map((calcProc, index) => (
                <TableRow key={index}>
                  <TableCell>
                    <Typography>
                      {calcProc.mainfile}
                    </Typography>
                    <Typography variant="caption">
                      {calcProc.calc_hash}
                    </Typography>
                  </TableCell>
                  <TableCell>
                    <Typography>
                      {calcProc.parser_name}
                    </Typography>
                  </TableCell>
                  <TableCell>
                    <Typography>
                      {calcProc.current_task_name}
                    </Typography>
                    <Typography variant="caption">
                      task&nbsp;
                      <b>
                        [{calcProc.task_names.indexOf(calcProc.current_task_name) + 1}/{calcProc.task_names.length}]
                      </b>
                    </Typography>
                  </TableCell>
                </TableRow>
              ))}
            </Table>
          )}
          <div className={classes.detailsContent}>
            <ReactJson src={upload} enableClipboard={false} collapsed={1} />
          </div>
        </ExpansionPanelDetails>
      </ExpansionPanel>
    )
  }
}

export default withStyles(Upload.styles)(Upload);