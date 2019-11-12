import React from 'react'
import Button from '@material-ui/core/Button'
import TextField from '@material-ui/core/TextField'
import Dialog from '@material-ui/core/Dialog'
import DialogActions from '@material-ui/core/DialogActions'
import DialogContent from '@material-ui/core/DialogContent'
import DialogContentText from '@material-ui/core/DialogContentText'
import DialogTitle from '@material-ui/core/DialogTitle'
import PropTypes from 'prop-types'
import { IconButton, Tooltip, withStyles, Paper, MenuItem, Popper, CircularProgress } from '@material-ui/core'
import EditIcon from '@material-ui/icons/Edit'
import AddIcon from '@material-ui/icons/Add'
import RemoveIcon from '@material-ui/icons/Delete'
import Autosuggest from 'react-autosuggest'
import match from 'autosuggest-highlight/match'
import parse from 'autosuggest-highlight/parse'
import { compose } from 'recompose'
import { withApi } from './api'

class SuggestionsTextFieldUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    value: PropTypes.string.isRequired,
    onChange: PropTypes.func.isRequired,
    suggestions: PropTypes.func.isRequired,
    suggestionValue: PropTypes.func.isRequired,
    suggestionRendered: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {
      width: '100%'
    },
    suggestion: {
      display: 'block'
    },
    suggestionsList: {
      margin: 0,
      padding: 0,
      listStyleType: 'none'
    },
    popper: {
      zIndex: theme.zIndex.modal + 200
    }
  })

  constructor(props) {
    super(props)
    this.lastRequestId = null
    this.unmounted = false
  }

  componentWillUnmount() {
    this.unmounted = true
  }

  componentDidMount() {
    this.unmounted = false
  }

  loadSuggestions(value) {
    if (this.state.isLoading) {
      return
    }

    if (this.lastRequestId !== null) {
      clearTimeout(this.lastRequestId)
    }

    this.setState({
      isLoading: true
    })

    this.lastRequestId = setTimeout(() => {
      this.props.suggestions(value).then(suggestions => {
        if (!this.unmounted) {
          this.setState({
            isLoading: false,
            suggestions: suggestions
          })
        }
      })
    }, 200)
  }

  state = {
    suggestions: [],
    isLoading: false,
    anchorEl: null
  }

  renderInputComponent(inputProps) {
    const { classes, inputRef = () => {}, ref, ...other } = inputProps

    return (
      <TextField
        fullWidth
        InputProps={{
          inputRef: node => {
            ref(node)
            inputRef(node)
          },
          classes: {
            input: classes.input
          }
        }}
        {...other}
      />
    )
  }

  renderSuggestion(suggestion, { query, isHighlighted }) {
    suggestion = this.props.suggestionRendered(suggestion)
    const matches = match(suggestion, query)
    const parts = parse(suggestion, matches)

    return (
      <MenuItem selected={isHighlighted} component="div">
        <div>
          {parts.map(part => (
            <span key={part.text} style={{ fontWeight: part.highlight ? 500 : 400 }}>
              {part.text}
            </span>
          ))}
        </div>
      </MenuItem>
    )
  }

  render() {
    const { classes, onChange, value, suggestions, suggestionValue, suggestionRendered, ...props } = this.props
    const { anchorEl } = this.state

    const handleSuggestionsFetchRequested = ({ value }) => {
      this.loadSuggestions(value)
    }

    const handleSuggestionsClearRequested = () => {
      this.setState({suggestions: []})
    }

    const handleChange = (event, { newValue }) => {
      onChange({target: {value: newValue}})
    }

    const autosuggestProps = {
      renderInputComponent: this.renderInputComponent.bind(this),
      suggestions: this.state.suggestions,
      onSuggestionsFetchRequested: handleSuggestionsFetchRequested,
      onSuggestionsClearRequested: handleSuggestionsClearRequested,
      getSuggestionValue: suggestionValue,
      renderSuggestion: this.renderSuggestion.bind(this)
    }

    return (
      <div className={classes.root}>
        <Autosuggest
          {...autosuggestProps}
          inputProps={{
            classes,
            value: value,
            onChange: handleChange,
            inputRef: node => {
              this.setState({anchorEl: node})
            },
            InputLabelProps: {
              shrink: true
            },
            ...props
          }}
          theme={{
            suggestionsList: classes.suggestionsList,
            suggestion: classes.suggestion
          }}
          renderSuggestionsContainer={options => (
            <Popper anchorEl={anchorEl} open={Boolean(options.children)} className={classes.popper}>
              <Paper
                square
                {...options.containerProps}
                style={{ width: anchorEl ? anchorEl.clientWidth : undefined }}
              >
                {options.children}
              </Paper>
            </Popper>
          )}
        />
      </div>
    )
  }
}

const SuggestionsTextField = withStyles(SuggestionsTextFieldUnstyled.styles)(SuggestionsTextFieldUnstyled)

var urlPattern = new RegExp('^(https?:\\/\\/)?' + // protocol
  '((([a-z\\d]([a-z\\d-]*[a-z\\d])*)\\.?)+[a-z]{2,}|' + // domain name
  '((\\d{1,3}\\.){3}\\d{1,3}))' + // OR ip (v4) address
  '(\\:\\d+)?(\\/[-a-z\\d%_.~+]*)*' + // port and path
  '(\\?[;&a-z\\d%_.~+=-]*)?' + // query string
  '(\\#[-a-z\\d_]*)?$', 'i') // fragment locator

function isURL(str) {
  return str === '' || urlPattern.test(str.trim())
}

class ListTextInputUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    values: PropTypes.arrayOf(PropTypes.object).isRequired,
    validate: PropTypes.func,
    label: PropTypes.string,
    errorLabel: PropTypes.string,
    onChange: PropTypes.func,
    component: PropTypes.func
  }

  static styles = theme => ({
    root: {},
    row: {
      display: 'flex'
    },
    buttonContainer: {
      position: 'relative',
      width: 52
    },
    button: {
      position: 'absolute',
      bottom: 0
    }
  })

  render() {
    const { classes, values, onChange, label, errorLabel, validate, component, ...fieldProps } = this.props
    const handleChange = (index, value) => {
      // TODO
      if (onChange) {
        const newValues = [...values]
        if (newValues[index]) {
          newValues[index].value = value
        } else {
          newValues[index] = {value: value}
        }
        onChange(newValues)
      }
    }

    const handleAdd = () => {
      if (onChange) {
        onChange([...values, {value: ''}])
      }
    }

    const handleRemove = (index) => {
      if (onChange) {
        onChange([...values.slice(0, index), ...values.slice(index + 1)])
      }
    }

    const Component = component || TextField
    const normalizedValues = values.length === 0 ? [{value: ''}] : values

    return <React.Fragment>
      {normalizedValues.map(({value, message, success}, index) => {
        let error = validate && !validate(value)
        let labelValue
        if (index === 0) {
          labelValue = label
        }
        if (error) {
          labelValue = errorLabel || 'Bad value'
        } else if (message) {
          labelValue = message
          error = !success
        }
        return <div key={index} className={classes.row}>
          <Component
            value={value}
            error={error}
            onChange={(event) => handleChange(index, event.target.value)}
            label={labelValue}
            margin={index === 0 ? 'normal' : 'dense'}
            InputLabelProps={{
              shrink: true
            }}
            {...fieldProps}
          />
          <div className={classes.buttonContainer}>
            {normalizedValues.length > 1
              ? <IconButton className={classes.button} size="tiny" onClick={() => handleRemove(index)}>
                <RemoveIcon fontSize="inherit" />
              </IconButton> : ''}
          </div>
          <div className={classes.buttonContainer}>
            {index + 1 === normalizedValues.length && normalizedValues[index].value !== ''
              ? <IconButton className={classes.button} size="tiny" onClick={handleAdd}>
                <AddIcon fontSize="inherit" />
              </IconButton> : ''}
          </div>
        </div>
      })}
    </React.Fragment>
  }
}

const ListTextInput = withStyles(ListTextInputUnstyled.styles)(ListTextInputUnstyled)

class SuggestionsListTextInput extends React.Component {
  render() {
    return <ListTextInput component={SuggestionsTextField} {...this.props} />
  }
}

class EditUserMetadataDialogUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    query: PropTypes.object.isRequired,
    total: PropTypes.number,
    example: PropTypes.object,
    buttonProps: PropTypes.object,
    api: PropTypes.object.isRequired,
    raiseError: PropTypes.func.isRequired,
    user: PropTypes.object,
    onEditComplete: PropTypes.func,
    disabled: PropTypes.bool,
    title: PropTypes.string
  }

  static styles = theme => ({
    dialog: {
      width: '100%'
    },
    submitWrapper: {
      margin: theme.spacing.unit,
      position: 'relative'
    },
    submitProgress: {
      position: 'absolute',
      top: '50%',
      left: '50%',
      marginTop: -12,
      marginLeft: -12
    }
  })

  constructor(props) {
    super(props)
    this.handleButtonClick = this.handleButtonClick.bind(this)
    this.handleClose = this.handleClose.bind(this)
    this.handleSubmit = this.handleSubmit.bind(this)
    this.verifyTimer = null
    this.state = {...this.defaultState}
    this.editData = {
      comment: '',
      references: [],
      coauthors: [],
      shared_with: [],
      datasets: [],
      with_embargo: true
    }
    this.unmounted = false
  }

  defaultState = {
    open: false,
    actions: {},
    isVerifying: false,
    verified: true,
    submitting: false
  }

  componentWillUnmount() {
    this.unmounted = true
  }

  update() {
    const { example } = this.props
    this.editData = {
      comment: example.comment || '',
      references: example.references || [],
      coauthors: (example.authors || []).filter(author => author.user_id !== example.uploader.user_id).map(author => author.email),
      shared_with: (example.owners || []).filter(author => author.user_id !== example.uploader.user_id).map(author => author.email),
      datasets: (example.datasets || []).map(ds => ds.name),
      with_embargo: example.with_embargo
    }
  }

  componentDidMount() {
    this.unmounted = false
    this.update()
  }

  componentDidUpdate(prevProps) {
    if (prevProps.example.calc_id !== this.props.example.calc_id) {
      this.update()
    }
  }

  verify() {
    if (this.state.isVerifying) {
      return
    }

    if (this.verifyTimer !== null) {
      clearTimeout(this.verifyTimer)
    }

    this.setState({
      isVerifying: true, verified: false
    })

    this.verifyTimer = setTimeout(() => {
      this.submitPromise(true).then(newState => {
        this.setState(newState)
      }).catch(error => {
        this.setState({verified: false, isVerifying: false})
        return this.props.raiseError(error)
      })
    }, 200)
  }

  submitPromise(verify) {
    const { query, api } = this.props
    const { actions } = this.state

    const editRequest = {
      query: query,
      verify: verify,
      actions: actions
    }
    return api.edit(editRequest).then(data => {
      if (this.unmounted) {
        return
      }

      const newActions = {...this.state.actions}
      let verified = true
      if (data.actions) {
        Object.keys(newActions).forEach(key => {
          if (Array.isArray(newActions[key])) {
            newActions[key] = newActions[key].map((action, i) => {
              verified &= !data.actions[key] || data.actions[key].success !== false
              return data.actions[key]
                ? {...(data.actions[key][i] || {}), value: action.value}
                : action
            })
          }
        })
      }
      return {actions: newActions, isVerifying: false, verified: verified}
    })
  }

  handleButtonClick() {
    const { open } = this.state
    if (!open) {
      this.update()
    }

    this.setState({open: !open})
  }

  handleClose() {
    this.setState({submitting: true})
    this.setState({...this.defaultState})
  }

  handleSubmit() {
    this.setState({submitting: true})

    this.submitPromise(false).then(newState => {
      if (this.props.onEditComplete) {
        this.props.onEditComplete()
      }
      this.setState({...newState, submitting: false})
      this.handleClose()
    }).catch(error => {
      this.setState({verified: false, isVerifying: false, submitting: false})
      return this.props.raiseError(error)
    })
  }

  render() {
    const { classes, buttonProps, total, api, user, example, disabled, title } = this.props
    const { open, actions, verified, submitting } = this.state

    const dialogEnabled = user && example.uploader && example.uploader.user_id === user.sub && !disabled
    const submitEnabled = Object.keys(actions).length && !submitting && verified

    const listTextInputProps = (key, verify) => {
      const values = actions[key] ? actions[key] : this.editData[key].map(value => ({value: value}))
      return {
        id: key,
        fullWidth: true,
        values: values,
        onChange: values => {
          this.setState({actions: {...actions, [key]: values}})
          if (verify) {
            this.verify()
          }
        }
      }
    }

    const userSuggestions = query => {
      return api.getUsers(query)
        .then(result => result.users)
        .catch(err => {
          console.error(err)
          return []
        })
    }

    return (
      <React.Fragment>
        <IconButton {...(buttonProps || {})} onClick={this.handleButtonClick} disabled={!dialogEnabled}>
          <Tooltip title={title || `Edit user metadata${dialogEnabled ? '' : '. You can only edit your data.'}`}>
            <EditIcon />
          </Tooltip>
        </IconButton>
        {dialogEnabled
          ? <Dialog classes={{paper: classes.dialog}} open={open} onClose={this.handleClose} disableBackdropClick disableEscapeKeyDown>
            <DialogTitle>Edit the user metadata of {total} entries</DialogTitle>
            <DialogContent>
              <DialogContentText>
                You are editing {total} {total === 1 ? 'entry' : 'entries'}. {total > 1
                  ? 'The fields are pre-filled with data from the first entry for.' : ''
                } Only the fields that you change will be updated.
                Be aware that all references, co-authors, shared_with, or datasets count as
                one field.
              </DialogContentText>
              <TextField
                // id="comment"
                label="Comment"
                value={actions.comment !== undefined ? actions.comment.value : this.editData.comment}
                onChange={event => this.setState({actions: {...actions, comment: {value: event.target.value}}})}
                margin="normal"
                multiline rows="4"
                fullWidth
                placeholder="Add a comment"
                InputLabelProps={{ shrink: true }}
              />
              <ListTextInput
                label="References"
                {...listTextInputProps('references')}
                errorLabel="References must be valid URLs"
                placeholder="Add a URL reference"
                validate={isURL}
              />
              <SuggestionsListTextInput
                label="Co-authors"
                suggestions={userSuggestions}
                suggestionValue={v => {
                  return v.email
                }}
                suggestionRendered={v => `${v.name} (${v.email})`}
                placeholder="Add a co-author by name"
                {...listTextInputProps('coauthors', true)}
              />
              <SuggestionsListTextInput
                {...listTextInputProps('shared_with', true)}
                suggestions={userSuggestions}
                suggestionValue={v => v.email}
                suggestionRendered={v => `${v.name} (${v.email})`}
                label="Shared with"
                placeholder="Add a user by name to share with"
              />
              <SuggestionsListTextInput
                {...listTextInputProps('datasets', true)}
                suggestions={prefix => {
                  return api.getDatasets(prefix)
                    .then(result => result.results.map(ds => ds.name))
                    .catch(err => {
                      console.error(err)
                      return []
                    })
                }}
                suggestionValue={v => v}
                suggestionRendered={v => v}
                label="Datasets"
                placeholder="Add a dataset"
              />
            </DialogContent>
            {Object.keys(actions).length
              ? <DialogContent>
                <DialogContentText>
                    The following fields will be updated with the given values: <i>
                    {Object.keys(actions).map(action => action).join(', ')}</i>.
                    Updating many entries might take a few seconds.
                </DialogContentText>
              </DialogContent>
              : ''}
            <DialogActions>
              <Button onClick={this.handleClose} color="primary" disabled={submitting}>
                Cancel
              </Button>
              <div className={classes.submitWrapper}>
                <Button onClick={this.handleSubmit} color="primary" disabled={!submitEnabled}>
                  Submit
                </Button>
                {submitting && <CircularProgress size={24} className={classes.submitProgress} />}
              </div>
            </DialogActions>
          </Dialog>
          : ''
        }
      </React.Fragment>
    )
  }
}

export default compose(withApi(false, false), withStyles(EditUserMetadataDialogUnstyled.styles))(EditUserMetadataDialogUnstyled)
