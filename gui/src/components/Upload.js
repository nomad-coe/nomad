import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, ExpansionPanel, ExpansionPanelSummary, Typography,
  ExpansionPanelDetails, Stepper, Step, StepLabel, Table, TableRow, TableCell, TableBody,
  Checkbox,
  FormControlLabel,
  TablePagination,
  TableHead,
  Tooltip,
  LinearProgress,
  CircularProgress} from '@material-ui/core'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ReactJson from 'react-json-view'
import CalcLinks from './CalcLinks'
import { compose } from 'recompose'
import { withErrors } from './errors'
import { debug } from '../config'

class Upload extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    upload: PropTypes.object.isRequired,
    checked: PropTypes.bool,
    onCheckboxChanged: PropTypes.func
  }

  static styles = theme => ({
    root: {},
    heading: {
      fontSize: theme.typography.pxToRem(15),
      fontWeight: theme.typography.fontWeightRegular
    },
    details: {
      padding: 0,
      display: 'block',
      overflowX: 'auto'
    },
    summary: {
      overflowX: 'auto'
    },
    detailsContent: {
      margin: theme.spacing.unit * 3
    },
    title: {
      flexBasis: '20%',
      flexShrink: 0,
      marginRight: theme.spacing.unit * 2
    },
    checkbox: {
      marginRight: theme.spacing.unit * 2
    },
    stepper: {
      width: '100%',
      padding: 0
    },
    buttonCell: {
      overflow: 'hidden',
      whiteSpace: 'nowrap',
      textAlign: 'right'
    },
    progress: {
      marginLeft: -theme.spacing.unit * 0.5,
      width: theme.spacing.unit * 13 - 2,
      alignItems: 'center',
      display: 'flex'
    }
  });

  state = {
    upload: this.props.upload,
    page: 1,
    rowsPerPage: 5
  }

  updateUpload() {
    window.setTimeout(() => {
      this.state.upload.update()
        .then(upload => {
          console.assert(upload.proc, 'Uploads always must have a proc')
          this.setState({upload: upload})
          if (upload.proc.status !== 'SUCCESS' && upload.proc.status !== 'FAILURE' && !upload.proc.is_stale) {
            this.updateUpload()
          }
        })
        .catch(error => {
          this.setState({upload: null})
          this.props.raiseError(error)
        })
    }, 500)
  }

  componentDidMount() {
    this.updateUpload()
  }

  handleChangePage = (_, page) => {
    this.setState({page: page + 1})
  }

  handleChangeRowsPerPage = event => {
    const rowsPerPage = event.target.value
    this.setState({rowsPerPage: rowsPerPage})
  }

  onCheckboxChanged(_, checked) {
    if (this.props.onCheckboxChanged) {
      this.props.onCheckboxChanged(checked)
    }
  }

  renderTitle() {
    const { classes } = this.props
    const { name, upload_id, create_time } = this.state.upload

    return (
      <div className={classes.title}>
        <Typography variant="title">
          {name || upload_id}
        </Typography>
        <Typography variant="subheading">
          {new Date(Date.parse(create_time)).toLocaleString()}
        </Typography>
      </div>
    )
  }

  renderStepper() {
    const { classes } = this.props
    const { calc_procs, task_names, current_task_name, status, errors } = this.state.upload.proc

    let activeStep = task_names.indexOf(current_task_name)
    activeStep += (status === 'SUCCESS') ? 1 : 0

    const labelPropsFactories = {
      parse_all: (props) => {
        props.children = 'parse'
        if (calc_procs.length > 0) {
          const failures = calc_procs.filter(calcProc => calcProc.status === 'FAILURE')
          if (failures.length) {
            props.error = true
            props.optional = (
              <Typography variant="caption" color="error">
                {calc_procs.filter(p => p.status === 'SUCCESS').length}/{calc_procs.length}
                , {failures.length} failed
              </Typography>
            )
          } else {
            props.optional = (
              <Typography variant="caption">
                {calc_procs.filter(p => p.status === 'SUCCESS').length}/{calc_procs.length}
              </Typography>
            )
          }

        } else if (status === 'SUCCESS') {
          props.error = true
          props.optional = (
            <Typography variant="caption" color="error">No calculations found.</Typography>
          )
        }
      }
    }

    return (
      <Stepper activeStep={activeStep} classes={{root: classes.stepper}}>
        {task_names.map((label, index) => {
          const labelProps = {
            children: label,
            error: activeStep === index && status === 'FAILURE'
          }

          const labelPropsFactory = labelPropsFactories[label]
          if (labelPropsFactory) {
            labelPropsFactory(labelProps)
          }

          if (labelProps.error && status === 'FAILURE') {
            labelProps.optional = (
              <Typography variant="caption" color="error">
                {errors.join(' ')}
              </Typography>
            )
          }

          return (
            <Step key={label}>
              <StepLabel {...labelProps} />
            </Step>
          )
        })}
      </Stepper>
    )
  }

  renderCalcTable() {
    const { classes } = this.props
    const { page, rowsPerPage } = this.state
    const { calc_procs, status, upload_hash } = this.state.upload.proc

    if (calc_procs.length === 0) {
      if (this.state.upload.is_ready) {
        return (
          <Typography className={classes.detailsContent}>
            {status === 'SUCCESS' ? 'No calculcations found.' : 'There are errors and no calculations to show.'}
          </Typography>
        )
      } else {
        return (
          <Typography className={classes.detailsContent}>
            Processing ...
          </Typography>
        )
      }
    }

    const renderRow = (calcProc, index) => {
      const { mainfile, calc_hash, parser_name, task_names, current_task_name, status, errors } = calcProc
      const color = status === 'FAILURE' ? 'error' : 'default'
      const row = (
        <TableRow key={index}>
          <TableCell>
            <Typography color={color}>
              {mainfile}
            </Typography>
            <Typography variant="caption" color={color}>
              {calc_hash}
            </Typography>
          </TableCell>
          <TableCell>
            <Typography color={color}>
              {parser_name.replace('parsers/', '')}
            </Typography>
          </TableCell>
          <TableCell>
            <Typography color={color}>
              {current_task_name}
            </Typography>
            <Typography variant="caption" color={color}>
              task&nbsp;
              <b>
                [{task_names.indexOf(current_task_name) + 1}/{task_names.length}]
              </b>
            </Typography>
          </TableCell>
          <TableCell>
            {status === 'SUCCESS'
              ? <CalcLinks uploadHash={upload_hash} calcHash={calc_hash} /> : ''}
          </TableCell>
        </TableRow>
      )

      if (status === 'FAILURE') {
        return (
          <Tooltip title={errors.map((error, index) => (<p key={index}>{error}</p>))}>
            {row}
          </Tooltip>
        )
      } else {
        return row
      }
    }

    const total = calc_procs.length
    const emptyRows = rowsPerPage - Math.min(rowsPerPage, total - (page - 1) * rowsPerPage)

    return (
      <Table>
        <TableHead>
          <TableRow>
            <TableCell>mainfile</TableCell>
            <TableCell>code</TableCell>
            <TableCell>task</TableCell>
            <TableCell></TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          {calc_procs.slice((page - 1) * rowsPerPage, page * rowsPerPage).map(renderRow)}
          {emptyRows > 0 && (
            <TableRow style={{ height: 57 * emptyRows }}>
              <TableCell colSpan={6} />
            </TableRow>
          )}
          <TableRow>
            <TablePagination
              count={total}
              rowsPerPage={rowsPerPage}
              page={page - 1}
              onChangePage={this.handleChangePage}
              onChangeRowsPerPage={this.handleChangeRowsPerPage}
            />
          </TableRow>
        </TableBody>
      </Table>
    )
  }

  render() {
    const { classes } = this.props
    const { upload } = this.state

    if (this.state.upload) {
      return (
        <ExpansionPanel>
          <ExpansionPanelSummary
            expandIcon={<ExpandMoreIcon/>} classes={{root: classes.summary}}>
            {!upload.is_ready
              ?
                <div className={classes.progress}>
                  <CircularProgress size={32}/>
                </div>
              :
                <FormControlLabel control={(
                  <Checkbox
                    checked={this.props.checked}
                    className={classes.checkbox}
                    onClickCapture={(e) => e.stopPropagation()}
                    onChange={this.onCheckboxChanged.bind(this)}
                  />
                )}/>
            }
            {this.renderTitle()} {this.renderStepper()}
          </ExpansionPanelSummary>
          <ExpansionPanelDetails style={{width: '100%'}} classes={{root: classes.details}}>
            {this.renderCalcTable()}
            {debug
              ? <div className={classes.detailsContent}>
                <ReactJson src={upload} enableClipboard={false} collapsed={0} />
              </div> : ''}
          </ExpansionPanelDetails>
        </ExpansionPanel>
      )
    } else {
      return ''
    }
  }
}

export default compose(withErrors, withStyles(Upload.styles))(Upload)
