import React from 'react'
import PropTypes from 'prop-types'
import { withStyles, ExpansionPanel, ExpansionPanelSummary, Typography,
  ExpansionPanelDetails, Stepper, Step, StepLabel,
  Checkbox, FormControlLabel, Tooltip,
  CircularProgress} from '@material-ui/core'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ReactJson from 'react-json-view'
import { compose } from 'recompose'
import { withErrors } from '../errors'
import { withRouter } from 'react-router'
import { debug } from '../../config'
import EntryList, { EntryListUnstyled } from '../search/EntryList'
import { withDomain } from '../domains'

class Upload extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    upload: PropTypes.object.isRequired,
    checked: PropTypes.bool,
    onCheckboxChanged: PropTypes.func,
    onDoesNotExist: PropTypes.func,
    onPublished: PropTypes.func,
    history: PropTypes.any.isRequired,
    domain: PropTypes.object.isRequired,
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
    shortTitle: {
      textOverflow: 'ellipsis',
      whiteSpace: 'nowrap',
      overflowX: 'inherit',
      direction: 'rtl',
      textAlign: 'left'
    },
    title: {
      textOverflow: 'ellipsis',
      whiteSpace: 'nowrap',
      overflowX: 'inherit'
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
  })

  static defaultSelectedColumns = ['mainfile', 'parser', 'proc', 'tasks_status']

  state = {
    upload: this.props.upload,
    params: {
      page: 1,
      per_page: 10,
      order_by: 'tasks_status',
      order: 1
    },
    updating: true, // it is still not complete and continuously looking for updates
    columns: {}
  }

  _unmounted = false

  constructor(props) {
    super(props)
    this.handleChange = this.handleChange.bind(this)
  }

  componentDidUpdate(prevProps, prevState) {
    if (prevProps.domain != this.props.domain) {
      this.updateColumns()
    }

    if (this.state.updating) {
      return
    }

    if (this.state.params === prevState.params && prevProps.upload.process_running === this.props.upload.process_running) {
      return
    }

    this.update()
  }

  updateColumns() {
    const { domain } = this.props

    const domainColumns = domain ? domain.searchResultColumns : {}
    const otherColumns = {...domainColumns, ...EntryListUnstyled.defaultColumns}
    Object.keys(otherColumns).forEach(key => {
      otherColumns[key] = {
        ...otherColumns[key],
        supportsSort: false
      }
    })
    const columns = {
      mainfile: {
        label: 'Mainfile',
        supportsSort: true,
        description: 'The path to the main output of this entry in the upload.'
      },
      parser: {
        label: 'Parser',
        supportsSort: true,
        description: 'The parser that was used to process this entry.',
        render: entry => entry.parser.replace('parsers/', '')
      },
      proc: {
        label: 'Processing',
        supportsSort: false,
        description: 'Details on the processing of this entry.',
        render: entry => `${entry.current_task} [${entry.tasks.indexOf(entry.current_task) + 1}/${entry.tasks.length}]`
      },
      tasks_status: {
        label: 'Status',
        supportsSort: true,
        descriptions: 'Entry processing status',
        render: entry => {
          const { tasks_status, errors, warnings } = entry
          const label = tasks_status.toLowerCase()
          const error = tasks_status === 'FAILURE' || errors.length > 0 || warnings.length > 0
          let tooltip = null
          if (tasks_status === 'FAILURE') {
            tooltip = `Calculation processing failed with errors: ${errors.join(', ')}`
          }
          if (errors.length > 0) {
            tooltip = `Calculation processed with errors: ${errors.join(', ')}`
          }
          if (warnings.length > 0) {
            tooltip = `Calculation processed with warnings: ${warnings.join(', ')}`
          }

          if (error) {
            return <Tooltip title={tooltip}>
                <Typography color="error">
                  {label}
                </Typography>
              </Tooltip>
          } else {
            return label
          }
        }
      },
      ...otherColumns
    }
    this.setState({columns: columns})
  }

  update() {
    if (this._unmounted) {
      return
    }

    const {page, per_page, order_by, order} = this.state.params
    const wasPublished = this.state.published
    this.state.upload.get(page, per_page, order_by, order)
      .then(upload => {
        const {tasks_running, process_running, current_task, published} = upload
        if (!this._unmounted) {
          if (published && !wasPublished) {
            if (this.props.onPublished) {
              this.props.onPublished()
            }
          }
          const continueUpdating = tasks_running || process_running || current_task === 'uploading'
          this.setState({upload: upload, updating: continueUpdating})
          if (continueUpdating) {
            window.setTimeout(() => {
              this.update()
            }, 500)
          }
        }
      })
      .catch(error => {
        if (!this._unmounted) {
          if (error.name === 'DoesNotExist') {
            this.props.onDoesNotExist()
          } else {
            this.props.raiseError(error)
          }
        }
      })
  }

  componentDidMount() {
    this.updateColumns()
    this.update()
  }

  componentWillUnmount() {
    this._unmounted = true
  }

  handleChange(changes) {
    this.setState({params: {...this.state.params, ...changes}})
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
        <Typography variant="h6" className={name ? classes.shortTitle : classes.title}>
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
    const { columns } = this.state
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

    const data = {
      pagination: calcs.pagination,
      results: calcs.results.map(calc => ({
        ...calc.metadata, ...calc
      }))
    }

    return <EntryList
      query={{upload_id: this.props.upload_id}}
      columns={columns}
      selectedColumns={Upload.defaultSelectedColumns}
      editable
      data={data}
      onChange={this.handleChange}
      {...this.state.params}
    />
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

export default compose(withRouter, withErrors, withDomain, withStyles(Upload.styles))(Upload)
