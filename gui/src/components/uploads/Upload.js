import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, ExpansionPanel, ExpansionPanelSummary, Typography,
  ExpansionPanelDetails, Stepper, Step, StepLabel, Table, TableRow, TableCell, TableBody,
  Checkbox, FormControlLabel, TablePagination, TableHead, Tooltip,
  CircularProgress,
  TableSortLabel} from '@material-ui/core'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ReactJson from 'react-json-view'
import { compose } from 'recompose'
import { withErrors } from '../errors'
import { withRouter } from 'react-router'
import { debug } from '../../config'

class Upload extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    upload: PropTypes.object.isRequired,
    checked: PropTypes.bool,
    onCheckboxChanged: PropTypes.func,
    onDoesNotExist: PropTypes.func,
    onPublished: PropTypes.func,
    history: PropTypes.any.isRequired
  }

  static styles = theme => ({
    root: {
      marginBottom: theme.spacing.unit
    },
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
    titleContainer: {
      flex: '0 0 auto',
      marginRight: theme.spacing.unit * 2,
      width: 350,
      overflowX: 'hidden'
    },
    title: {
      textOverflow: 'ellipsis',
      whiteSpace: 'nowrap',
      overflowX: 'inherit',
      direction: 'rtl',
      textAlign: 'left'
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
    clickableRow: {
      cursor: 'pointer'
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
    updating: true // it is still not complete and continuously looking for updates
  }

  _unmounted = false

  update(params) {
    if (this._unmounted) {
      return
    }

    const {page, perPage, orderBy, order} = params
    this.state.upload.get(page, perPage, orderBy, order === 'asc' ? 1 : -1)
      .then(upload => {
        const {tasks_running, process_running, current_task, published} = upload
        if (!this._unmounted) {
          if (published) {
            this.setState({...params})
            if (this.props.onPublished) {
              this.props.onPublished()
            }
            return
          }
          const continueUpdating = tasks_running || process_running || current_task === 'uploading'
          this.setState({upload: upload, params: params, updating: continueUpdating})
          if (continueUpdating) {
            window.setTimeout(() => {
              this.update(this.state.params)
            }, 500)
          }
        }
      })
      .catch(error => {
        if (!this._unmounted) {
          this.setState({...params})
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
      <div className={classes.titleContainer}>
        <Typography variant="h6" className={classes.title}>
          {name || new Date(Date.parse(create_time)).toLocaleString()}
        </Typography>
        {name
          ? <Typography variant="subtitle1">
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
    const { calcs, tasks, current_task, tasks_running, tasks_status, process_running, current_process } = upload

    // map tasks [ uploading, extracting, parse_all, cleanup ] to steps
    const steps = [ 'upload', 'process', 'publish' ]
    let step = null
    const task_index = tasks.indexOf(current_task)
    if (task_index === 0) {
      step = 'upload'
    } else if (task_index > 0 && tasks_running) {
      step = 'process'
    } else if (!upload.published) {
      step = 'publish'
    }
    const stepIndex = upload.published ? steps.length : steps.indexOf(step)

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
              processing failed
            </Typography>
          )
        }
      },
      publish: (props) => {
        if (upload.published) {
          props.children = 'published'
        } else {
          props.children = 'inspect'

          if (process_running) {
            if (current_process === 'publish_upload') {
              props.children = 'approved'
              props.optional = <Typography variant="caption">moving data ...</Typography>
            } else if (current_process === 'delete_upload') {
              props.children = 'declined'
              props.optional = <Typography variant="caption">deleting data ...</Typography>
            }
          } else {
            props.optional = <Typography variant="caption">publish or delete</Typography>
          }
        }
      }
    }

    return (
      <Stepper activeStep={stepIndex} classes={{root: classes.stepper}}>
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
      const { mainfile, upload_id, calc_id, parser, tasks, current_task, tasks_status, errors } = calc
      let tooltip_html = null
      let color = tasks_status === 'FAILURE' ? 'error' : 'default'
      if (tasks_status === 'FAILURE') {
        tooltip_html = `Calculation processing failed with errors: ${errors.join(', ')}`
        color = 'error'
      }
      if (calc.errors.length > 0) {
        color = 'error'
        tooltip_html = `Calculation processed with errors: ${calc.errors.join(', ')}`
      }
      if (calc.warnings.length > 0) {
        color = 'error'
        tooltip_html = `Calculation processed with warnings: ${calc.warnings.join(', ')}`
      }
      const processed = tasks_status === 'FAILURE' || tasks_status === 'SUCCESS'
      const row = (
        <TableRow key={index} hover={processed}
          onClick={() => this.props.history.push(`/uploads/${upload_id}/${calc_id}`)}
          className={processed ? classes.clickableRow : null} >

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
            <Typography color={color}>{tasks_status.toLowerCase()}</Typography>
          </TableCell>
        </TableRow>
      )

      if (tooltip_html) {
        return (
          <Tooltip key={calc_id} title={tooltip_html}>
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
      { id: 'tasks_status', sort: true, label: 'status' }
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

  renderCheckBox() {
    const { classes } = this.props
    const { upload } = this.state

    if (upload.tasks_running || upload.process_running) {
      return <div className={classes.progress}>
        <CircularProgress size={32}/>
      </div>
    } else if (!upload.published) {
      return <FormControlLabel control={(
        <Checkbox
          checked={this.props.checked}
          className={classes.checkbox}
          onClickCapture={(e) => e.stopPropagation()}
          onChange={this.onCheckboxChanged.bind(this)}
        />
      )}/>
    } else {
      return ''
    }
  }

  render() {
    const { classes } = this.props
    const { upload } = this.state
    const { errors } = upload

    if (this.state.upload) {
      return (
        <div className={classes.root}>
          <ExpansionPanel>
            <ExpansionPanelSummary
              expandIcon={<ExpandMoreIcon/>}
              classes={{root: classes.summary}}>

              {this.renderCheckBox()} {this.renderTitle()} {this.renderStepper()}
            </ExpansionPanelSummary>
            <ExpansionPanelDetails style={{width: '100%'}} classes={{root: classes.details}}>
              {errors && errors.length > 0
                ? <Typography className={classes.detailsContent} color="error">
                  Upload processing has errors: {errors.join(', ')}
                </Typography> : ''
              }
              {upload.calcs ? this.renderCalcTable() : ''}
              {debug
                ? <div className={classes.detailsContent}>
                  <ReactJson src={upload} enableClipboard={false} collapsed={0} />
                </div> : ''}
            </ExpansionPanelDetails>
          </ExpansionPanel>
        </div>
      )
    } else {
      return ''
    }
  }
}

export default compose(withRouter, withErrors, withStyles(Upload.styles))(Upload)
