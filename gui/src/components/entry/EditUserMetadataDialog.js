/*
 * Copyright The NOMAD Authors.
 *
 * This file is part of NOMAD. See https://nomad-lab.eu for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
import React from 'react'
import Button from '@material-ui/core/Button'
import TextField from '@material-ui/core/TextField'
import Dialog from '@material-ui/core/Dialog'
import DialogActions from '@material-ui/core/DialogActions'
import DialogContent from '@material-ui/core/DialogContent'
import DialogContentText from '@material-ui/core/DialogContentText'
import DialogTitle from '@material-ui/core/DialogTitle'
import PropTypes from 'prop-types'
import { IconButton, Tooltip, withStyles, Paper, MenuItem, Popper, CircularProgress,
  FormGroup, Checkbox } from '@material-ui/core'
import EditIcon from '@material-ui/icons/Edit'
import AddIcon from '@material-ui/icons/Add'
import RemoveIcon from '@material-ui/icons/Delete'
import Autosuggest from 'react-autosuggest'
import match from 'autosuggest-highlight/match'
import parse from 'autosuggest-highlight/parse'
import { withApi } from '../api'

const local_users = {}

function update_local_user(user) {
  local_users[user.user_id] = user
}

class MyAutosuggestUnstyled extends React.PureComponent {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    shouldRenderSuggestions: PropTypes.func,
    suggestions: PropTypes.func.isRequired,
    getSuggestionValue: PropTypes.func.isRequired,
    getSuggestionRenderValue: PropTypes.func,
    inputProps: PropTypes.object,
    value: PropTypes.any,
    onChange: PropTypes.func,
    allowNew: PropTypes.bool
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

  state = {
    suggestions: [],
    loadingRequest: null,
    inputValue: this.props.getSuggestionValue(this.props.value)
  }

  componentDidUpdate(prevProps) {
    if (prevProps.value !== this.props.value && this.props.value !== undefined) {
      this.setState({inputValue: this.props.getSuggestionValue(this.props.value)})
    }
  }

  popperNode = null
  unmounted = false
  lastRequested = null
  lastRequestTimeout = null

  componentWillUnmount() {
    this.unmounted = true
  }

  componentDidMount() {
    this.unmounted = false
  }

  handleSuggestionsFetchRequested({value}) {
    value = value.trim()
    this.lastRequested = value

    if (this.state.loadingRequest) {
      return
    }

    if (this.lastRequestTimeout !== null) {
      clearTimeout(this.lastRequestTimeout)
    }

    this.lastRequestTimeout = setTimeout(() => {
      this.setState({
        loadingRequest: this.lastRequested
      })
      this.props.suggestions(value).then(suggestions => {
        if (!this.unmounted) {
          if (this.lastRequested !== this.state.loadingRequest) {
            this.handleSuggestionsFetchRequested({value: this.lastRequested})
          }
          this.setState({
            loadingRequest: null,
            suggestions: suggestions
          }, () => this.handleChangedInputValue(this.state.inputValue))
        }
      })
    }, 200)
  }

  handleSuggestionsClearRequested() {
    this.setState({suggestions: []})
  }

  handleChangedInputValue(value) {
    const {allowNew} = this.props
    const normalizedValue = value.trim().toLowerCase()
    const {getSuggestionValue, onChange} = this.props
    const event = value => ({target: {value: value}})
    if (onChange) {
      if (normalizedValue.length === 0) {
        onChange(event(null))
      } else {
        const {suggestions} = this.state
        const matchingSuggestion = suggestions
          .find(suggestion => getSuggestionValue(suggestion).toLowerCase() === normalizedValue)
        if (matchingSuggestion) {
          const matchingSuggestionValue = getSuggestionValue(matchingSuggestion)
          if (getSuggestionValue(matchingSuggestionValue) !== normalizedValue) {
            this.setState({inputValue: matchingSuggestionValue + value.trimLeft().slice(matchingSuggestionValue.length)})
          }
          onChange(event(matchingSuggestion))
        } else {
          if (allowNew) {
            onChange(event(value.trim()))
          } else {
            onChange(event(undefined))
          }
        }
      }
    }
  }

  handleChange(event, { newValue }) {
    this.setState({inputValue: newValue})
    this.handleChangedInputValue(newValue)
  }

  renderSuggestion(suggestion, { query, isHighlighted }) {
    const getValue = this.props.getSuggestionRenderValue || this.props.getSuggestionValue
    const inputValue = getValue(suggestion)
    const matches = match(inputValue, query)
    const parts = parse(inputValue, matches)
    return (
      <MenuItem selected={isHighlighted} component="div">
        <div>
          {parts.map((part, i) => (
            <span key={i} style={{ fontWeight: part.highlight ? 500 : 300 }}>
              {part.text}
            </span>
          ))}
        </div>
      </MenuItem>
    )
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
          name: 'search', // try to prevent browsers ignore autocomplete="off"
          type: 'search', // try to prevent browsers ignore autocomplete="off"
          classes: {
            input: classes.input
          }
        }}
        {...other}
      />
    )
  }

  renderSuggestionsContainer(options) {
    const {classes} = this.props
    return <Popper anchorEl={this.popperNode} open={Boolean(options.children)} className={classes.popper}>
      <Paper
        square
        {...options.containerProps}
        style={{ width: this.popperNode ? this.popperNode.clientWidth : null }}
      >
        {options.children}
      </Paper>
    </Popper>
  }

  render() {
    const {
      classes, shouldRenderSuggestions, getSuggestionValue, suggestions, allowNew,
      getSuggestionRenderValue, value, onChange, ...inputProps} = this.props

    return <div className={classes.root}>
      <Autosuggest
        renderInputComponent={this.renderInputComponent.bind(this)}
        renderSuggestion={this.renderSuggestion.bind(this)}
        suggestions={this.state.suggestions}
        onSuggestionsFetchRequested={this.handleSuggestionsFetchRequested.bind(this)}
        onSuggestionsClearRequested={this.handleSuggestionsClearRequested.bind(this)}
        getSuggestionValue={getSuggestionValue}
        shouldRenderSuggestions={shouldRenderSuggestions}
        inputProps={{
          classes,
          value: this.state.inputValue,
          onChange: this.handleChange.bind(this),
          inputRef: node => {
            this.popperNode = node
          },
          InputLabelProps: {
            shrink: true
          },
          ...inputProps
        }}
        theme={{
          suggestionsList: classes.suggestionsList,
          suggestion: classes.suggestion
        }}
        renderSuggestionsContainer={this.renderSuggestionsContainer.bind(this)}
      />
    </div>
  }
}

const MyAutosuggest = withStyles(MyAutosuggestUnstyled.styles)(MyAutosuggestUnstyled)

class DatasetInputUnstyled extends React.Component {
  static propTypes = {
    value: PropTypes.string, // name
    label: PropTypes.string,
    error: PropTypes.bool,
    api: PropTypes.object.isRequired,
    onChange: PropTypes.func,
    margin: PropTypes.any
  }

  suggestions(query) {
    const {api} = this.props
    query = query.toLowerCase()
    return api.getDatasets(query)
      .then(result => result.results.map(ds => ds.dataset_name))
      .catch(err => {
        console.error(err)
        return []
      })
  }

  getSuggestionRenderValue(suggestion) {
    return suggestion
  }

  getSuggestionValue(suggestion) {
    return suggestion || ''
  }

  render() {
    const {label, onChange, value, margin} = this.props
    let usedLabel = label
    if (value === undefined) {
      usedLabel = 'This dataset does not exist, it will be created'
    }

    return <MyAutosuggest onChange={onChange} value={value}
      allowNew
      suggestions={this.suggestions.bind(this)}
      getSuggestionValue={this.getSuggestionValue.bind(this)}
      getSuggestionRenderValue={this.getSuggestionRenderValue.bind(this)}
      shouldRenderSuggestions={() => true}
      margin={margin}
      label={usedLabel}
      placeholder="Type the dataset's name"
    />
  }
}

const DatasetInput = withApi(DatasetInputUnstyled)

class ReferenceInput extends React.Component {
  static propTypes = {
    onChange: PropTypes.func.isRequired,
    value: PropTypes.string,
    label: PropTypes.string
  }

  state = {
    inputValue: this.props.value || ''
  }

  componentDidUpdate(prevProps) {
    if (prevProps.value !== this.props.value && this.props.value !== undefined) {
      this.setState({inputValue: this.props.value || ''})
    }
  }

  handleChange(event) {
    const inputValue = event.target.value
    this.setState({inputValue: inputValue})
    const trimmedInputValue = inputValue.trim()
    let value = null
    if (trimmedInputValue.length !== 0) {
      if (isURL(trimmedInputValue)) {
        value = trimmedInputValue
      } else {
        value = undefined
      }
    }
    if (value !== this.props.value) {
      this.props.onChange({target: {value: value}})
    }
  }

  render() {
    const {value, onChange, label, ...rest} = this.props
    return <TextField
      fullWidth
      {...rest}
      type="search" name="search" // attempt to avoid browsers autofill, since they seem to ignore autocomplete="off"
      value={this.state.inputValue}
      onChange={this.handleChange.bind(this)}
      error={value === undefined}
      label={value === undefined ? 'A reference must be a valid url' : label}
      placeholder="Enter a URL to a related resource"
    />
  }
}

class ActionInput extends React.PureComponent {
  static propTypes = {
    value: PropTypes.object.isRequired,
    onChange: PropTypes.func.isRequired,
    label: PropTypes.string,
    component: PropTypes.elementType
  }

  handleChange(event) {
    const value = event.target.value
    if (value !== this.props.value.value) {
      this.props.onChange({massage: null, success: true, value: value})
    }
  }

  render() {
    const {value, onChange, component, label, ...rest} = this.props
    const Component = component || TextField

    const {message, success} = value
    let labelWithMessageAndError = label
    let error = false
    if (success === false) {
      labelWithMessageAndError = message || 'Bad value'
      error = true
    } else if (message) {
      labelWithMessageAndError = message
    }

    return <Component
      value={value.value}
      label={labelWithMessageAndError}
      error={error}
      onChange={this.handleChange.bind(this)}
      {...rest} />
  }
}

var urlPattern = new RegExp('^(https?:\\/\\/)?' + // protocol
  '((([a-z\\d]([a-z\\d-]*[a-z\\d])*)\\.?)+[a-z]{2,}|' + // domain name
  '((\\d{1,3}\\.){3}\\d{1,3}))' + // OR ip (v4) address
  '(\\:\\d+)?(\\/[-a-z\\d%_.~+]*)*' + // port and path
  '(\\?[;&a-z\\d%_.~+=-]*)?' + // query string
  '(\\#[-a-z\\d_]*)?$', 'i') // fragment locator

function isURL(str) {
  return !str || str === '' || urlPattern.test(str.trim())
}

class ListTextInputUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    values: PropTypes.arrayOf(PropTypes.object).isRequired,
    label: PropTypes.string,
    onChange: PropTypes.func,
    component: PropTypes.any
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
    const { classes, values, onChange, label, component, ...fieldProps } = this.props

    const handleChange = (index, value) => {
      if (onChange) {
        const newValues = [...values]
        newValues[index] = value
        onChange(newValues)
      }
    }

    const handleAdd = () => {
      if (onChange) {
        onChange([...values, {value: null}])
      }
    }

    const handleRemove = (index) => {
      if (onChange) {
        onChange([...values.slice(0, index), ...values.slice(index + 1)])
      }
    }

    const Component = component || TextField
    const normalizedValues = values.length === 0 ? [{value: null}] : values

    return <React.Fragment>
      {normalizedValues.map((value, index) => {
        let labelValue
        if (index === 0) {
          labelValue = label
        }
        return <div key={index} className={classes.row}>
          <ActionInput component={Component}
            value={value}
            onChange={value => handleChange(index, value)}
            label={labelValue}
            margin={index === 0 ? 'normal' : 'dense'}
            InputLabelProps={{
              shrink: true
            }}
            {...fieldProps}
          />
          <div className={classes.buttonContainer}>
            {normalizedValues.length > 1 || (normalizedValues.length === 1 && normalizedValues[0].value)
              ? <IconButton className={classes.button} size="medium" onClick={() => handleRemove(index)}>
                <RemoveIcon fontSize="inherit" />
              </IconButton> : ''}
          </div>
          <div className={classes.buttonContainer}>
            {index + 1 === normalizedValues.length && value.value
              ? <IconButton className={classes.button} size="medium" onClick={handleAdd}>
                <AddIcon fontSize="inherit" />
              </IconButton> : ''}
          </div>
        </div>
      })}
    </React.Fragment>
  }
}

const ListTextInput = withStyles(ListTextInputUnstyled.styles)(ListTextInputUnstyled)

class InviteUserDialogUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    api: PropTypes.object.isRequired
  }

  static styles = theme => ({
    button: {
      marginLeft: theme.spacing(1)
    },
    dialog: {
      width: '100%'
    },
    submitWrapper: {
      margin: theme.spacing(1),
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

  defaultState = {
    open: false,
    data: {
      first_name: '',
      last_name: '',
      email: '',
      affiliation: ''
    },
    error: null,
    submitting: false,
    submitEnabled: false
  }

  state = {...this.defaultState}

  handleClose() {
    this.setState({...this.defaultState, open: false})
  }

  handleSubmit() {
    this.setState({submitting: true})

    this.props.api.inviteUser(this.state.data).then(() => {
      this.handleClose()
    }).catch(error => {
      const detail = error?.response?.data?.detail
      if (detail) {
        this.setState({error: detail, submitting: false, submitEnabled: false})
      } else {
        this.setState({error: '' + error, submitting: false, submitEnabled: false})
      }
    })
  }

  handleChange(key, value) {
    const {data} = this.state
    const valid = value && !Object.keys(data).find(dataKey => !(key === dataKey || data[dataKey]))
    this.setState({submitEnabled: valid, data: {...data, [key]: value}})
  }

  handleOpen() {
    this.setState({
      ...this.defaultState, open: true
    })
  }

  render() {
    const {classes} = this.props
    const {open, data, submitting, submitEnabled, error} = this.state
    const input = (key, label) => <TextField
      label={label}
      value={data[key]}
      onChange={event => this.handleChange(key, event.target.value)}
      margin="normal"
      fullWidth
    />
    return <React.Fragment>
      <Button className={classes.button}
        onClick={this.handleOpen.bind(this)}
        color="secondary" disabled={submitting}
      >
        Invite new user
      </Button>
      <Dialog
        classes={{paper: classes.dialog}}
        open={open}
        onClose={this.handleClose.bind(this)} disableBackdropClick disableEscapeKeyDown>
        <DialogTitle>Invite a new user to NOMAD</DialogTitle>
        <DialogContent>
          <DialogContentText>
            If you want to add a user as co-author or share your data with someone that
            is not already a NOMAD user, you can invite this person here. We need just a few
            details about this person. After your invite, the new user will receive an
            Email that allows her to set a password and further details. Anyhow, you will
            be able to add the user as co-author or someone to share with immediately after the
            invite.
          </DialogContentText>
          {error && <DialogContentText color="error">
            {error}
          </DialogContentText>}
          {input('email', 'Email')}
          {input('first_name', 'First name')}
          {input('last_name', 'Last name')}
          {input('affiliation', 'Affiliation')}
        </DialogContent>
        <DialogActions>
          <Button onClick={this.handleClose.bind(this)} disabled={submitting}>
            Cancel
          </Button>
          <div className={classes.submitWrapper}>
            <Button onClick={this.handleSubmit.bind(this)} color="primary" disabled={!submitEnabled}>
              Submit
            </Button>
            {submitting && <CircularProgress size={24} className={classes.submitProgress} />}
          </div>
        </DialogActions>
      </Dialog>
    </React.Fragment>
  }
}

const InviteUserDialog = withApi(withStyles(InviteUserDialogUnstyled.styles)(InviteUserDialogUnstyled))

class UserMetadataFieldUnstyled extends React.PureComponent {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    children: PropTypes.node,
    modified: PropTypes.bool,
    onChange: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {
      flexWrap: 'nowrap',
      alignItems: 'flex-start',
      marginTop: theme.spacing(2)
    },
    container: {
      width: '100%'
    },
    checkbox: {
      marginLeft: -theme.spacing(1) * 2,
      marginRight: theme.spacing(1),
      marginTop: theme.spacing(1)
    }
  })

  render() {
    const {children, classes, modified, onChange} = this.props
    return <FormGroup row className={classes.root}>
      <Checkbox
        classes={{root: classes.checkbox}}
        checked={modified}
        onChange={(event, checked) => onChange(checked)}
      />
      <div className={classes.container}>
        {children}
      </div>
    </FormGroup>
  }
}

const UserMetadataField = withStyles(UserMetadataFieldUnstyled.styles)(UserMetadataFieldUnstyled)

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
    title: PropTypes.string,
    info: PropTypes.object,
    text: PropTypes.string
  }

  static styles = theme => ({
    dialog: {
      width: '100%'
    },
    error: {
      marginTop: theme.spacing(2)
    },
    submitWrapper: {
      margin: theme.spacing(1),
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
      entry_coauthors: [],
      reviewers: [],
      datasets: []
    }
    this.unmounted = false
  }

  defaultState = {
    open: false,
    actions: {},
    success: true,
    message: null,
    isVerifying: false,
    verified: true,
    submitting: false
  }

  componentWillUnmount() {
    this.unmounted = true
  }

  update() {
    const { example } = this.props
    if (example.authors) {
      example.authors.forEach(user => update_local_user(user))
    }
    if (example.viewers) {
      example.viewers.forEach(user => update_local_user(user))
    }
    this.editData = {
      comment: example.comment || '',
      references: example.references || [],
      entry_coauthors: (example.authors || [])
        .filter(user => user.user_id !== example.main_author.user_id)
        .map(user => user.user_id),
      reviewers: (example.viewers || [])
        .filter(user => user.user_id !== example.main_author.user_id)
        .map(user => user.user_id),
      datasets: (example.datasets || []).map(ds => ds.dataset_name)
    }
  }

  componentDidMount() {
    this.unmounted = false
    this.update()
  }

  componentDidUpdate(prevProps) {
    if (prevProps.example.entry_id !== this.props.example.entry_id) {
      this.update()
    }
  }

  verify() {
    if (this.state.isVerifying) {
      return
    }

    const { actions } = this.state

    for (let key of Object.keys(actions)) {
      const actionValues = actions[key]
      if (Array.isArray(actionValues)) {
        if (actionValues.find(action => action.value === undefined)) {
          this.setState({
            isVerifying: false, verified: false
          })
          return
        }
      }
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

    // remove null values to allow swagger validation
    const actionsCopy = {...actions}
    Object.keys(actionsCopy).forEach(key => {
      if (Array.isArray(actionsCopy[key])) {
        actionsCopy[key] = actionsCopy[key].map(action => {
          const actionCopy = {...action}
          if (!actionCopy.value) {
            delete actionCopy.value
          }
          return actionCopy
        })
      }
    })

    const editRequest = {
      owner: query.owner,
      query: { ...query, owner: undefined },
      verify: verify,
      actions: actionsCopy
    }

    return api.edit(editRequest).then(data => {
      if (this.unmounted) {
        return
      }

      let verified = true
      const newActions = {...this.state.actions}
      if (data.actions) {
        Object.keys(newActions).forEach(key => {
          if (Array.isArray(newActions[key])) {
            newActions[key] = newActions[key].map((action, i) => {
              verified &= !data.actions[key] || data.actions[key][i].success !== false
              return data.actions[key]
                ? {...(data.actions[key][i] || {}), value: action.value}
                : action
            })
          }
        })
      }
      return {
        actions: newActions,
        success: data.success,
        message: data.message,
        isVerifying: false,
        verified: verified && data.success
      }
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
    const { classes, buttonProps, total, user, example, disabled, title } = this.props
    const { open, actions, verified, submitting, success, message } = this.state

    const dialogEnabled = user && example.main_author && example.main_author.user_id === user.sub && !disabled
    const submitEnabled = Object.keys(actions).length && !submitting && verified

    const editDataToActions = editData => {
      if (Array.isArray(editData)) {
        return editData.map(value => ({value: value}))
      } else {
        return {value: editData}
      }
    }

    const listTextInputProps = (key, verify) => {
      const values = actions[key] ? actions[key] : editDataToActions(this.editData[key])

      return {
        values: values,
        onChange: values => {
          this.setState({actions: {...actions, [key]: values}}, () => {
            if (verify) {
              this.verify()
            }
          })
        }
      }
    }

    const metadataFieldProps = (key, verify, defaultValue) => ({
      modified: Boolean(actions[key]),
      defaultValue: defaultValue,
      onChange: checked => {
        if (checked) {
          this.setState({actions: {...actions, [key]: editDataToActions(this.editData[key] || defaultValue)}}, () => {
            if (verify) {
              this.verify()
            }
          })
        } else {
          this.setState({actions: {...actions, [key]: undefined}})
        }
      }
    })

    const tooltipProps = {
      title: title || `Edit user metadata${dialogEnabled ? '' : '. You can only edit your data.'}`
    }
    const allButtonProps = {
      onClick: this.handleButtonClick,
      disabled: !dialogEnabled,
      ...(buttonProps || {})
    }

    return (
      <React.Fragment>
        {!this.props.text &&
          <IconButton {...allButtonProps}>
            <Tooltip {...tooltipProps}>
              <EditIcon />
            </Tooltip>
          </IconButton>}
        {this.props.text &&
          <Tooltip {...tooltipProps}>
            <Button {...allButtonProps}>
              {this.props.text}
            </Button>
          </Tooltip>}
        {dialogEnabled
          ? <Dialog classes={{paper: classes.dialog}} open={open} onClose={this.handleClose} disableBackdropClick disableEscapeKeyDown>
            <DialogTitle>Edit the user metadata of {total} entries</DialogTitle>
            <DialogContent>
              <DialogContentText>
                You are editing {total} {total === 1 ? 'entry' : 'entries'}. {total > 1
                  ? 'The fields are pre-filled with data from the first entry for.' : ''
                } Only the checked fields will be updated.
                The fields references, co-authors, shared with users,
                and datasets can have many values. Changing one value, will apply all values.
              </DialogContentText>
              {!success && <DialogContentText color="error" className={classes.error}>
                {message}
              </DialogContentText>}
              <UserMetadataField {...metadataFieldProps('comment')}>
                <ActionInput component={TextField}
                  label="Comment"
                  value={actions.comment !== undefined ? actions.comment : {value: this.editData.comment}}
                  onChange={value => this.setState({actions: {...actions, comment: value}})}
                  margin="normal"
                  multiline
                  rowsMax="10"
                  fullWidth
                  placeholder="Add a comment"
                  InputLabelProps={{ shrink: true }}
                />
              </UserMetadataField>
              <UserMetadataField {...metadataFieldProps('references', true)}>
                <ListTextInput
                  component={ReferenceInput}
                  {...listTextInputProps('references', true)}
                  label="References"
                />
              </UserMetadataField>
              <UserMetadataField {...metadataFieldProps('datasets', true)}>
                <ListTextInput
                  component={DatasetInput}
                  {...listTextInputProps('datasets', true)}
                  label="Datasets"
                />
              </UserMetadataField>
            </DialogContent>
            {this.renderDialogActions(submitting, submitEnabled)}
          </Dialog>
          : ''
        }
      </React.Fragment>
    )
  }

  renderDialogActions(submitting, submitEnabled) {
    const {classes, info} = this.props

    if (submitting) {
      return <DialogActions>
        <DialogContentText style={{marginLeft: 16}}>Do not close the page. This might take up to several minutes for editing many entries.</DialogContentText>
        <span style={{flexGrow: 1}} />
        <div className={classes.submitWrapper}>
          <Button onClick={this.handleSubmit} disabled={!submitEnabled} color="primary">
            Submit
          </Button>
          <CircularProgress size={24} className={classes.submitProgress} />
        </div>
      </DialogActions>
    } else {
      return <DialogActions>
        {info && !info.oasis && <InviteUserDialog />}
        <span style={{flexGrow: 1}} />
        <Button onClick={this.handleClose} disabled={submitting}>
          Cancel
        </Button>
        <Button onClick={this.handleSubmit} disabled={!submitEnabled} color="primary">
          Submit
        </Button>
      </DialogActions>
    }
  }
}

export default withApi(withStyles(EditUserMetadataDialogUnstyled.styles)(EditUserMetadataDialogUnstyled))
