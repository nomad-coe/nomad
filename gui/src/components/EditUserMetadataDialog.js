import React from 'react'
import Button from '@material-ui/core/Button'
import TextField from '@material-ui/core/TextField'
import Dialog from '@material-ui/core/Dialog'
import DialogActions from '@material-ui/core/DialogActions'
import DialogContent from '@material-ui/core/DialogContent'
import DialogContentText from '@material-ui/core/DialogContentText'
import DialogTitle from '@material-ui/core/DialogTitle'
import PropTypes from 'prop-types'
import { IconButton, Tooltip, withStyles, Paper, MenuItem, Popper } from '@material-ui/core'
import EditIcon from '@material-ui/icons/Edit'
import AddIcon from '@material-ui/icons/Add'
import RemoveIcon from '@material-ui/icons/Delete'
import ReactJson from 'react-json-view'
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
        this.setState({
          isLoading: false,
          suggestions: suggestions
        })
      })
    }, 1000)
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
    values: PropTypes.arrayOf(PropTypes.string).isRequired,
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
        onChange([...values.slice(0, index), value, ...values.slice(index + 1)])
      }
    }

    const handleAdd = () => {
      if (onChange) {
        onChange([...values, ''])
      }
    }

    const handleRemove = (index) => {
      if (onChange) {
        onChange([...values.slice(0, index), ...values.slice(index + 1)])
      }
    }

    const Component = component || TextField
    const normalizedValues = values.length === 0 ? [''] : values

    return <React.Fragment>
      {normalizedValues.map((value, index) => {
        const error = validate && !validate(value)
        let labelValue = index === 0 ? label : null
        if (error) {
          labelValue = errorLabel || 'Bad value'
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
            {index + 1 === normalizedValues.length && normalizedValues[index] !== ''
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
    total: PropTypes.number,
    example: PropTypes.object,
    buttonProps: PropTypes.object,
    api: PropTypes.object.isRequired
  }

  static styles = theme => ({
    dialog: {
      width: '100%'
    }
  })

  constructor(props) {
    super(props)
    this.handleButtonClick = this.handleButtonClick.bind(this)
  }

  state = {
    open: false,
    editData: {
      comment: '',
      references: [],
      coAuthors: [],
      sharedWith: [],
      datasets: [],
      withEmbargo: true
    }
  }

  update() {
    const { example } = this.props
    const editData = {
      comment: example.comment || '',
      references: example.references || [],
      coAuthors: example.authors.filter(author => author.user_id !== example.uploader.user_id).map(author => author.email),
      sharedWith: example.owners.filter(author => author.user_id !== example.uploader.user_id).map(author => author.email),
      datasets: (example.datasets || []).map(ds => ds.name),
      withEmbargo: example.with_embargo
    }
    this.setState({editData: editData})
  }

  componentDidMount() {
    this.update()
  }

  componentDidUpdate(prevProps) {
    if (prevProps.example.calc_id !== this.props.example.calc_id) {
      this.update()
    }
  }

  handleButtonClick() {
    const { open } = this.state
    if (!open) {
      this.update()
    }

    this.setState({open: !open})
  }

  render() {
    const { classes, buttonProps, total, api } = this.props
    const { open } = this.state
    const close = () => this.setState({open: false})

    const handleChange = (key, value) => {
      this.setState({editData: {...this.state.editData, [key]: value}})
    }
    const value = key => this.state.editData[key]

    const userSuggestions = query => {
      return api.getUsers(query)
        .then(result => result.users)
        .catch((err) => {
          console.log(err)
          return []
        })
    }

    return (
      <React.Fragment>
        <Tooltip title="Edit user metadata">
          <IconButton {...(buttonProps || {})} onClick={this.handleButtonClick}>
            <EditIcon />
          </IconButton>
        </Tooltip>
        <Dialog classes={{paper: classes.dialog}} open={open} onClose={close} disableBackdropClick disableEscapeKeyDown>
          <DialogTitle>Edit the user metadata of {total} entries</DialogTitle>
          <DialogContent>
            <DialogContentText>
              TODO better text
            </DialogContentText>
            <TextField
              id="comment"
              label="Comment"
              value={value('comment')}
              onChange={event => handleChange('comment', event.target.value)}
              margin="normal"
              multiline
              fullWidth
            />
            <ListTextInput
              id="references"
              label="References"
              errorLabel="References must be valid URLs"
              placeholder="Add a URL reference"
              values={value('references')}
              onChange={values => handleChange('references', values)}
              validate={isURL}
              fullWidth
            />
            <SuggestionsListTextInput
              suggestions={userSuggestions}
              suggestionValue={v => v.email}
              suggestionRendered={v => `${v.name} (${v.email})`}
              id="coAuthors"
              label="Co-authors"
              placeholder="Add a co-author by name"
              values={value('coAuthors')}
              onChange={values => handleChange('coAuthors', values)}
              fullWidth
            />
            <SuggestionsListTextInput
              suggestions={userSuggestions}
              suggestionValue={v => v.email}
              suggestionRendered={v => `${v.name} (${v.email})`}
              id="sharedWith"
              label="Shared with"
              placeholder="Add a user by name to share with"
              values={value('sharedWith')}
              onChange={values => handleChange('sharedWith', values)}
              fullWidth
            />
            <SuggestionsListTextInput
              suggestions={prefix => {
                console.log(prefix)
                return api.getDatasets(prefix)
                  .then(result => result.results.map(ds => ds.name))
                  .catch((err) => {
                    console.log(err)
                    return []
                  })
              }}
              suggestionValue={v => v}
              suggestionRendered={v => v}
              id="datasets"
              label="Datasets"
              placeholder="Add a dataset"
              values={value('datasets')}
              onChange={values => handleChange('datasets', values)}
              fullWidth
            />
          </DialogContent>
          <DialogContent>
            <ReactJson
              src={this.state.editData}
              enableClipboard={false}
              collapsed={0}
            />
          </DialogContent>
          <DialogActions>
            <Button onClick={close} color="primary">
              Cancel
            </Button>
            <Button onClick={close} color="primary">
              Submit
            </Button>
          </DialogActions>
        </Dialog>
      </React.Fragment>
    )
  }
}

export default compose(withApi(false, false), withStyles(EditUserMetadataDialogUnstyled.styles))(EditUserMetadataDialogUnstyled)
