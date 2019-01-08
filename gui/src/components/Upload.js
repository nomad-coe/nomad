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
      orderBy: 'status',
      order: 'asc'
    },
    archiveLogs: null, // { uploadHash, calcHash } ids of archive to show logs for
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
        if (!this._unmounted) {
          const continueUpdating = upload.status !== 'SUCCESS' && upload.status !== 'FAILURE' && !upload.is_stale
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
          this.props.raiseError(error)
        }
      })
  }

  componentDidMount() {
    this.update(this.state.params)
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
    const { calcs, tasks, current_task, status, errors } = upload

    let activeStep = tasks.indexOf(current_task)
    activeStep += (status === 'SUCCESS') ? 1 : 0

    const labelPropsFactories = {
      uploading: (props) => {
        props.children = 'uploading'
        const { uploading } = upload
        if (upload.status !== 'FAILURE') {
          props.optional = (
            <Typography variant="caption">
              {uploading === 100  && current_task === tasks[0] ? 'waiting for processing ...' : `${uploading || 0}%`}
            </Typography>
          )
        }
      },
      extracting: (props) => {
        props.children = 'extracting'
        if (current_task === 'extracting') {
          props.optional = (
            <Typography variant="caption">
              be patient
            </Typography>
          )
        }
      },
      parse_all: (props) => {
        props.children = 'parse'
        if (!calcs) {
          props.optional = (
            <Typography variant="caption" >
              loading...
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
        {tasks.map((label, index) => {
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
    const { page, perPage, orderBy, order } = this.state.params
    const { calcs, status, waiting } = this.state.upload
    const { pagination, results } = calcs

    if (pagination.total === 0) {
      if (this.state.upload.completed) {
        return (
          <Typography className={classes.detailsContent}>
            {status === 'SUCCESS' ? 'No calculcations found.' : 'No calculations to show.'}
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
      const { mainfile, calc_hash, upload_hash, parser, tasks, current_task, status, errors } = calc
      const color = status === 'FAILURE' ? 'error' : 'default'
      const row = (
        <TableRow key={index}>
          <TableCell>
            <Typography color={color}>
              {mainfile}
            </Typography>
            <Typography variant="caption" color={color}>
              {upload_hash}/{calc_hash}
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
              {(status === 'SUCCESS' || status === 'FAILURE')
                ?
                  <a className={classes.logLink} href="#logs" onClick={() => this.setState({archiveLogs: { uploadHash: upload_hash, calcHash: calc_hash }})}>
                  {status.toLowerCase()}
                  </a>
                : status.toLowerCase()
              }
            </Typography>
          </TableCell>
          <TableCell>
            <CalcLinks uploadHash={upload_hash} calcHash={calc_hash} disabled={status !== 'SUCCESS'} />
          </TableCell>
        </TableRow>
      )

      if (status === 'FAILURE') {
        return (
          <Tooltip key={calc_hash} title={errors.map((error, index) => (<p key={`${calc_hash}-${index}`}>{error}</p>))}>
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
      { id: 'status', sort: true, label: 'status' },
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
          uploadHash={this.state.archiveLogs.uploadHash}
          calcHash={this.state.archiveLogs.calcHash}
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
              {!(upload.completed || upload.waiting)
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
