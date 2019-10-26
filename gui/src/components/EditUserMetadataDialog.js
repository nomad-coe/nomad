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
import deburr from 'lodash/deburr'

// TODO replace with the actual authors
const suggestions = [
  { label: 'Afghanistan' },
  { label: 'Aland Islands' },
  { label: 'Albania' },
  { label: 'Algeria' },
  { label: 'American Samoa' },
  { label: 'Andorra' },
  { label: 'Angola' },
  { label: 'Anguilla' },
  { label: 'Antarctica' },
  { label: 'Antigua and Barbuda' },
  { label: 'Argentina' },
  { label: 'Armenia' },
  { label: 'Aruba' },
  { label: 'Australia' },
  { label: 'Austria' },
  { label: 'Azerbaijan' },
  { label: 'Bahamas' },
  { label: 'Bahrain' },
  { label: 'Bangladesh' },
  { label: 'Barbados' },
  { label: 'Belarus' },
  { label: 'Belgium' },
  { label: 'Belize' },
  { label: 'Benin' },
  { label: 'Bermuda' },
  { label: 'Bhutan' },
  { label: 'Bolivia, Plurinational State of' },
  { label: 'Bonaire, Sint Eustatius and Saba' },
  { label: 'Bosnia and Herzegovina' },
  { label: 'Botswana' },
  { label: 'Bouvet Island' },
  { label: 'Brazil' },
  { label: 'British Indian Ocean Territory' },
  { label: 'Brunei Darussalam' }
]

class AuthorTextFieldUnstyled extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    value: PropTypes.string.isRequired,
    onChange: PropTypes.func.isRequired
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
    const matches = match(suggestion.label, query)
    const parts = parse(suggestion.label, matches)

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

  getSuggestions(value) {
    const inputValue = deburr(value.trim()).toLowerCase()
    const inputLength = inputValue.length
    let count = 0

    return inputLength === 0
      ? []
      : suggestions.filter(suggestion => {
        const keep =
            count < 5 && suggestion.label.slice(0, inputLength).toLowerCase() === inputValue

        if (keep) {
          count += 1
        }

        return keep
      })
  }

  getSuggestionValue(suggestion) {
    return suggestion.label
  }

  render() {
    const { classes, onChange, value, ...props } = this.props
    const { suggestions, anchorEl } = this.state

    const handleSuggestionsFetchRequested = ({ value }) => {
      this.setState({suggestions: this.getSuggestions(value)})
    }

    const handleSuggestionsClearRequested = () => {
      this.setState({suggestions: []})
    }

    const handleChange = (event, { newValue }) => {
      onChange({target: {value: newValue}})
    }

    const autosuggestProps = {
      renderInputComponent: this.renderInputComponent.bind(this),
      suggestions,
      onSuggestionsFetchRequested: handleSuggestionsFetchRequested,
      onSuggestionsClearRequested: handleSuggestionsClearRequested,
      getSuggestionValue: this.getSuggestionValue.bind(this),
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

const AuthorTextField = withStyles(AuthorTextFieldUnstyled.styles)(AuthorTextFieldUnstyled)

var urlPattern = new RegExp('^(https?:\\/\\/)?' + // protocol
  '((([a-z\\d]([a-z\\d-]*[a-z\\d])*)\\.?)+[a-z]{2,}|' + // domain name
  '((\\d{1,3}\\.){3}\\d{1,3}))' + // OR ip (v4) address
  '(\\:\\d+)?(\\/[-a-z\\d%_.~+]*)*' + // port and path
  '(\\?[;&a-z\\d%_.~+=-]*)?' + // query string
  '(\\#[-a-z\\d_]*)?$', 'i') // fragment locator

function isURL(str) {
  return urlPattern.test(str)
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

class AuthorsListTextInput extends React.Component {
  render() {
    return <ListTextInput component={AuthorTextField} {...this.props} />
  }
}

class EditUserMetadataDialog extends React.Component {
  static propTypes = {
    query: PropTypes.object
  }

  state = {
    open: false,
    comment: 'This is the existing comment and it is very long. Lorem ipsum dolor sit amet, consetetur sadipscing elitr, sed diam nonumy eirmod tempor invidunt ut labore et dolore magna aliquyam erat, sed diam voluptua. At vero eos et accusam et justo duo dolores et ea rebum. Stet clita kasd gubergren, no sea takimata sanctus est Lorem ipsum dolor sit amet.',
    references: ['http://reference1', 'http://reference2'],
    coAuthors: ['Scheidgen, Markus'],
    sharedWith: [],
    datasets: [],
    withEmbargo: true
  }

  handleChange(key, value) {
    this.setState({[key]: value})
  }

  render() {
    const { query, ...buttonProps } = this.props
    const { open } = this.state
    const close = () => this.setState({open: false})

    return (
      <React.Fragment>
        <Tooltip title="Edit user metadata">
          <IconButton {...buttonProps} onClick={() => this.setState({open: true})}>
            <EditIcon />
          </IconButton>
        </Tooltip>
        <Dialog open={open} onClose={close} disableBackdropClick disableEscapeKeyDown>
          <DialogTitle>Edit the user metadata of X entries</DialogTitle>
          <DialogContent>
            <DialogContentText>
              To subscribe to this website, please enter your email address here. We will send updates
              occasionally.
            </DialogContentText>
            <TextField
              id="comment"
              label="Comment"
              value={this.state.comment}
              onChange={event => this.handleChange('comment', event.target.value)}
              margin="normal"
              multiline
              fullWidth
            />
            <ListTextInput
              id="references"
              label="References"
              errorLabel="References must be valid URLs"
              placeholder="Add a URL reference"
              values={this.state.references}
              onChange={values => this.handleChange('references', values)}
              validate={isURL}
              fullWidth
            />
            <AuthorsListTextInput
              id="coAuthors"
              label="Co-authors"
              placeholder="Add a co-author by name"
              values={this.state.coAuthors}
              onChange={values => this.handleChange('coAuthors', values)}
              fullWidth
            />
            <AuthorsListTextInput
              id="sharedWith"
              label="Shared with"
              placeholder="Add a user by name to share with"
              values={this.state.sharedWith}
              onChange={values => this.handleChange('sharedWith', values)}
              fullWidth
            />
            <ListTextInput
              id="datasets"
              label="Datasets"
              placeholder="Add a dataset"
              values={this.state.datasets}
              onChange={values => this.handleChange('datasets', values)}
              fullWidth
            />
          </DialogContent>
          <DialogContent>
            <ReactJson src={this.props.query} enableClipboard={false} collapsed={0} />
          </DialogContent>
          <DialogActions>
            <Button onClick={close} color="primary">
              Cancel
            </Button>
            <Button onClick={close} color="primary">
              Subscribe
            </Button>
          </DialogActions>
        </Dialog>
      </React.Fragment>
    )
  }
}

export default EditUserMetadataDialog
