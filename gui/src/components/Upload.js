import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, ExpansionPanel, ExpansionPanelSummary, Typography,
  ExpansionPanelDetails, Stepper, Step, StepLabel, Table, TableRow, TableCell, TableBody,
  Checkbox, FormControlLabel, TablePagination, TableHead, Tooltip,
  CircularProgress,
  LinearProgress,
  TableSortLabel} from '@material-ui/core'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ReactJson from 'react-json-view'
import CalcLinks from './CalcLinks'
import { compose } from 'recompose'
import { withErrors } from './errors'
import { debug } from '../config'
import CalcProcLogPopper from './CalcProcLogPopper'
import CalcDialogButtons from './CalcDialogButtons'

class Upload extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    upload: PropTypes.object.isRequired,
    checked: PropTypes.bool,
    onCheckboxChanged: PropTypes.func,
    onDoesNotExist: PropTypes.func
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
    },
    logLink: {
      color: theme.palette.secondary.main,
      textDecoration: 'none',
      '&:hover': {
        textDecoration: 'underline'
      }
    }
  });

  state = {
    upload: this.props.upload,
    params: {
      page: 1,
      perPage: 5,
      orderBy: 'tasks_status',
      order: 'asc'
    },
    archiveLogs: null, // { uploadId, calcId } ids of archive to show logs for
    loading: true, // its loading data from the server and the user should know about it
    updating: true // it is still not complete and continieusly looking for updates
  }

  _unmounted = false

  update(params) {
    if (this._unmounted) {
      return
    }

    const {page, perPage, orderBy, order} = params
    this.setState({loading: true})
    this.state.upload.get(page, perPage, orderBy, order === 'asc' ? 1 : -1)
      .then(upload => {
        const {tasks_running, process_running, current_task} = upload
        if (!this._unmounted) {
          const continueUpdating = tasks_running || process_running || current_task === 'uploading'
          this.setState({upload: upload, loading: false, params: params, updating: continueUpdating})
          if (continueUpdating) {
            window.setTimeout(() => {
              if (!this.state.loading) {
                this.update(this.state.params)
              }
            }, 500)
          }
        }
      })
      .catch(error => {
        if (!this._unmounted) {
          this.setState({loading: false, ...params})
          if (error.name === 'DoesNotExist') {
            this.props.onDoesNotExist()
          } else {
            this.props.raiseError(error)
          }
        }
      })
  }

  componentDidMount() {
    this.update(this.state.params)
  }

  componentDidUpdate(prevProps) {
    if (!prevProps.upload.process_running && this.props.upload.process_running) {
      this.update(this.state.params)
    }
  }

  componentWillUnmount() {
    this._unmounted = true
  }

  handleChangePage = (_, page) => {
    this.update({...this.state.params, page: page + 1})
  }

  handleChangeRowsPerPage = event => {
    const perPage = event.target.value
    this.update({...this.state.params, perPage: perPage})
  }

  handleSort(orderBy) {
    let order = 'desc'
    if (this.state.params.orderBy === orderBy && this.state.params.order === 'desc') {
      order = 'asc'
    }
    this.update({...this.state.params, orderBy: orderBy, order: order})
  }

  onCheckboxChanged(_, checked) {
    if (this.props.onCheckboxChanged) {
      this.props.onCheckboxChanged(checked)
    }
  }

  renderTitle() {
    const { classes } = this.props
    const { name, create_time } = this.state.upload

    return (
      <div className={classes.title}>
        <Typography variant="title">
          {name || new Date(Date.parse(create_time)).toLocaleString()}
        </Typography>
        {name
          ? <Typography variant="subheading">
            {new Date(Date.parse(create_time)).toLocaleString()}
          </Typography>
          : 'this upload has no name'
        }
      </div>
    )
  }

  renderStepper() {
    const { classes } = this.props
    const { upload } = this.state
    const { calcs, tasks, current_task, tasks_running, tasks_status, process_running, current_process, errors } = upload

    // map tasks [ uploading, extracting, parse_all, cleanup ] to steps
    const steps = [ 'upload', 'process', 'commit' ]
    let step = null
    const task_index = tasks.indexOf(current_task)
    if (task_index === 0) {
      step = 'upload'
    } else if (task_index > 0 && tasks_running) {
      step = 'process'
    } else {
      step = 'commit'
    }
    const stepIndex = steps.indexOf(step)

    const labelPropsFactories = {
      upload: (props) => {
        if (step === 'upload') {
          props.children = 'uploading'
          const { uploading } = upload
          if (upload.tasks_status !== 'FAILURE') {
            props.optional = (
              <Typography variant="caption">
                {`${uploading || 0}%`}
              </Typography>
            )
          }
        } else {
          props.children = 'uploaded'
        }
      },
      process: (props) => {
        props.error = tasks_status === 'FAILURE'

        const processIndex = steps.indexOf('process')
        if (stepIndex <= processIndex) {
          props.children = 'processing'
        } else {
          props.children = 'processed'
        }

        if (current_task === 'extracting') {
          props.children = 'extracting'
          props.optional = (
            <Typography variant="caption">
              be patient
            </Typography>
          )
        } else if (current_task === 'parse_all') {
          props.children = 'parsing'
        }

        if (stepIndex >= processIndex) {
          if (!calcs) {
            props.optional = (
              <Typography variant="caption" >
                matching...
              </Typography>
            )
          } else if (calcs.pagination.total > 0) {
            const { total, successes, failures } = calcs.pagination
            if (failures) {
              props.error = true
              props.optional = (
                <Typography variant="caption" color="error">
                  {successes + failures}/{total}, {failures} failed
                </Typography>
              )
            } else {
              props.optional = (
                <Typography variant="caption">
                  {successes + failures}/{total}
                </Typography>
              )
            }
          } else if (tasks_status === 'SUCCESS') {
            props.error = true
            props.optional = (
              <Typography variant="caption" color="error">No calculations found.</Typography>
            )
          }
        }

        if (tasks_status === 'FAILURE') {
          props.optional = (
            <Typography variant="caption" color="error">
              {errors.join(' ')}
            </Typography>
          )
        }
      },
      commit: (props) => {
        props.children = 'inspect'

        if (process_running) {
          if (current_process === 'commit_upload') {
            props.children = 'approved'
            props.optional = <Typography variant="caption">moving data ...</Typography>
          } else if (current_process === 'delete_upload') {
            props.children = 'declined'
            props.optional = <Typography variant="caption">deleting data ...</Typography>
          }
        } else {
          props.optional = <Typography variant="caption">commit or delete</Typography>
        }
      }
    }

    return (
      <Stepper activeStep={steps.indexOf(step)} classes={{root: classes.stepper}}>
        {steps.map((label, index) => {
          const labelProps = {
            children: label
          }

          const labelPropsFactory = labelPropsFactories[label]
          if (labelPropsFactory) {
            labelPropsFactory(labelProps)
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
    const { page, perPage, orderBy, order } = this.state.params
    const { calcs, tasks_status, waiting } = this.state.upload
    const { pagination, results } = calcs

    if (pagination.total === 0) {
      if (!this.state.upload.tasks_running) {
        return (
          <Typography className={classes.detailsContent}>
            {tasks_status === 'SUCCESS' ? 'No calculcations found.' : 'No calculations to show.'}
          </Typography>
        )
      } else {
        if (waiting) {
          return (
            <Typography className={classes.detailsContent}>
                Uploading ...
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
    }

    const renderRow = (calc, index) => {
      const { mainfile, calc_id, upload_id, parser, tasks, current_task, tasks_status, errors } = calc
      const color = tasks_status === 'FAILURE' ? 'error' : 'default'
      const row = (
        <TableRow key={index}>
          <TableCell>
            <Typography color={color}>
              {mainfile}
            </Typography>
            <Typography variant="caption" color={color}>
              {upload_id}/{calc_id}
            </Typography>
          </TableCell>
          <TableCell>
            <Typography color={color}>
              {parser.replace('parsers/', '')}
            </Typography>
          </TableCell>
          <TableCell>
            <Typography color={color}>
              {current_task}
            </Typography>
            <Typography variant="caption" color={color}>
              task&nbsp;
              <b>
                [{tasks.indexOf(current_task) + 1}/{tasks.length}]
              </b>
            </Typography>
          </TableCell>
          <TableCell>
            <Typography color={color}>
              {(tasks_status === 'SUCCESS' || tasks_status === 'FAILURE')
                ? <a className={classes.logLink} href="#logs" onClick={() => this.setState({archiveLogs: { uploadId: upload_id, calcId: calc_id }})}>
                  {tasks_status.toLowerCase()}
                </a> : tasks_status.toLowerCase()
              }
            </Typography>
          </TableCell>
          <TableCell>
            <CalcDialogButtons uploadId={upload_id} calcId={calc_id} disabled={tasks_status !== 'SUCCESS'} />
          </TableCell>
        </TableRow>
      )

      if (tasks_status === 'FAILURE') {
        return (
          <Tooltip key={calc_id} title={errors.map((error, index) => (<p key={`${calc_id}-${index}`}>{error}</p>))}>
            {row}
          </Tooltip>
        )
      } else {
        return row
      }
    }

    const total = pagination.total
    const emptyRows = perPage - Math.min(perPage, total - (page - 1) * perPage)

    const columns = [
      { id: 'mainfile', sort: true, label: 'mainfile' },
      { id: 'parser', sort: true, label: 'code' },
      { id: 'task', sort: false, label: 'task' },
      { id: 'tasks_status', sort: true, label: 'status' },
      { id: 'links', sort: false, label: 'links' }
    ]

    return (
      <Table>
        <TableHead>
          <TableRow>
            {columns.map(column => (
              <TableCell key={column.id}>
                {column.sort
                  ? <Tooltip
                    title="Sort"
                    placement={'bottom-start'}
                    enterDelay={300}
                  >
                    <TableSortLabel
                      active={orderBy === column.id}
                      direction={order}
                      onClick={() => this.handleSort(column.id)}
                    >
                      {column.label}
                    </TableSortLabel>
                  </Tooltip>
                  : column.label
                }
              </TableCell>
            ))}
          </TableRow>
        </TableHead>
        <TableBody>
          {results.map(renderRow)}
          {emptyRows > 0 && (
            <TableRow style={{ height: 57 * emptyRows }}>
              <TableCell colSpan={6} />
            </TableRow>
          )}

          <TableRow>
            <TablePagination
              count={total}
              rowsPerPage={perPage}
              page={page - 1}
              onChangePage={this.handleChangePage}
              onChangeRowsPerPage={this.handleChangeRowsPerPage}
            />
          </TableRow>
        </TableBody>
      </Table>
    )
  }

  renderLogs() {
    if (this.state.archiveLogs) {
      return (
        <CalcProcLogPopper
          open={true}
          onClose={() => this.setState({archiveLogs: null})}
          anchorEl={window.parent.document.documentElement.firstElementChild}
          raiseError={this.props.raiseError}
          uploadId={this.state.archiveLogs.uploadId}
          calcId={this.state.archiveLogs.calcId}
        />
      )
    } else {
      return ''
    }
  }

  render() {
    const { classes } = this.props
    const { upload } = this.state

    if (this.state.upload) {
      return (
        <div ref={this.logPopperAnchor}>
          <ExpansionPanel>
            <ExpansionPanelSummary
              expandIcon={<ExpandMoreIcon/>} classes={{root: classes.summary}}>
              {(upload.tasks_running || upload.process_running)
                ? <div className={classes.progress}>
                  <CircularProgress size={32}/>
                </div>
                : <FormControlLabel control={(
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
              {upload.calcs ? this.renderCalcTable() : ''}
              {debug
                ? <div className={classes.detailsContent}>
                  <ReactJson src={upload} enableClipboard={false} collapsed={0} />
                </div> : ''}
              {this.state.loading && !this.state.updating ? <LinearProgress/> : ''}
            </ExpansionPanelDetails>
          </ExpansionPanel>
          {this.renderLogs()}
        </div>
      )
    } else {
      return ''
    }
  }
}

export default compose(withErrors, withStyles(Upload.styles))(Upload)
