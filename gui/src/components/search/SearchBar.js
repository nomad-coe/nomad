import React from 'react'
import PropTypes from 'prop-types'
import { withStyles } from '@material-ui/core/styles'
import ChipInput from 'material-ui-chip-input'
import Autosuggest from 'react-autosuggest'
import match from 'autosuggest-highlight/match'
import parse from 'autosuggest-highlight/parse'
import Paper from '@material-ui/core/Paper'
import MenuItem from '@material-ui/core/MenuItem'
import { Chip } from '@material-ui/core'
import { repoPrimaryColor } from '../../config'

function renderInput(inputProps) {
  const { classes, autoFocus, value, onChange, onAdd, onDelete, chips, ref, ...other } = inputProps

  return (
    <ChipInput
      clearInputValueOnChange
      onUpdateInput={onChange}
      onAdd={onAdd}
      onDelete={onDelete}
      value={chips}
      inputRef={ref}
      chipRenderer={
        ({ value, text, isFocused, isDisabled, handleClick, handleDelete, className }, key) => (
          <Chip
            key={key}
            className={className}
            style={{
              pointerEvents: isDisabled ? 'none' : undefined,
              backgroundColor: isFocused ? repoPrimaryColor[500] : undefined,
              color: isFocused ? 'white' : 'black'
            }}
            onClick={handleClick}
            onDelete={handleDelete}
            label={text}
          />
        )
      }
      {...other}
    />
  )
}

function renderSuggestion(suggestion, { query, isHighlighted }) {
  const matches = match(getSuggestionValue(suggestion), query)
  const parts = parse(getSuggestionValue(suggestion), matches)

  return (
    <MenuItem
      selected={isHighlighted}
      component='div'
      onMouseDown={(e) => e.preventDefault()} // prevent the click causing the input to be blurred
    >
      <div>
        {parts.map((part, index) => {
          return part.highlight ? (
            <span key={String(index)} style={{ fontWeight: 300 }}>
              {part.text}
            </span>
          ) : (
            <strong key={String(index)} style={{ fontWeight: 500 }}>
              {part.text}
            </strong>
          )
        })}
      </div>
    </MenuItem>
  )
}

function renderSuggestionsContainer(options) {
  const { containerProps, children } = options

  return (
    <Paper {...containerProps} square>
      {children}
    </Paper>
  )
}

function getSuggestionValue(suggestion) {
  return `${suggestion.key}=${suggestion.value}`
}

class SearchBar extends React.Component {
  static propTypes = {
    classes: PropTypes.object.isRequired,
    data: PropTypes.object.isRequired,
    searchValues: PropTypes.object.isRequired,
    onChanged: PropTypes.func.isRequired
  }

  static styles = theme => ({
    root: {},
    autosuggestRoot: {
      position: 'relative'
    },
    suggestionsContainerOpen: {
      position: 'absolute',
      zIndex: 100,
      marginTop: theme.spacing.unit,
      left: 0,
      right: 0
    },
    suggestion: {
      display: 'block'
    },
    suggestionsList: {
      margin: 0,
      padding: 0,
      listStyleType: 'none'
    },
    textField: {
      width: '100%'
    }
  })

  state = {
    suggestions: [],
    textFieldInput: ''
  }

  getSuggestions(value) {
    value = value.toLowerCase()

    const { data: { quantities } } = this.props
    const suggestions = []

    Object.keys(quantities).forEach(quantity => {
      Object.keys(quantities[quantity]).forEach(quantityValue => {
        if (quantityValue.toLowerCase().startsWith(value)) {
          suggestions.push({
            key: quantity,
            value: quantityValue
          })
        }
      })
    })

    return suggestions
  }

  handleSuggestionsFetchRequested = ({ value }) => {
    this.setState({
      suggestions: this.getSuggestions(value)
    })
  };

  handleSuggestionsClearRequested = () => {
    this.setState({
      suggestions: []
    })
  };

  handleTextFieldInputChange = (event, { newValue }) => {
    this.setState({
      textFieldInput: newValue
    })
  }

  handleAddChip(chip) {
    const values = {...this.props.searchValues}

    let key, value
    if (chip.includes('=')) {
      const parts = chip.split('=')
      key = parts[0]
      value = parts[1]
    } else {
      const suggestion = this.getSuggestions(chip)[0]
      key = suggestion.key
      value = suggestion.value
    }

    if (values[key]) {
      values[key] = key === 'atoms' ? [...values[key], value] : value
    } else {
      values[key] = key === 'atoms' ? [value] : value
    }

    this.setState({
      textFieldInput: ''
    })

    this.props.onChanged(values)
  }

  handleBeforeAddChip(chip) {
    const suggestions = this.getSuggestions(chip)
    if (suggestions.length > 0) {
      return true
    } else {
      return false
    }
  }

  handleDeleteChip(chip) {
    if (!chip) {
      return
    }
    const parts = chip.split('=')
    const key = parts[0]

    const values = {...this.props.searchValues}
    delete values[key]
    this.props.onChanged(values)
  }

  getChips() {
    const values = {...this.props.searchValues}
    return Object.keys(values).filter(key => values[key]).map(key => {
      if (key === 'atoms') {
        return `atoms=[${values[key].join(',')}]`
      } else {
        return `${key}=${values[key]}`
      }
    })
  }

  render() {
    const { classes, searchValues, onChanged, ...rest } = this.props

    return (
      <Autosuggest
        theme={{
          container: classes.autosuggestRoot,
          suggestionsContainerOpen: classes.suggestionsContainerOpen,
          suggestionsList: classes.suggestionsList,
          suggestion: classes.suggestion
        }}
        renderInputComponent={renderInput}
        suggestions={this.state.suggestions}
        onSuggestionsFetchRequested={this.handleSuggestionsFetchRequested}
        onSuggestionsClearRequested={this.handleSuggestionsClearRequested}
        renderSuggestionsContainer={renderSuggestionsContainer}
        getSuggestionValue={getSuggestionValue}
        renderSuggestion={renderSuggestion}
        onSuggestionSelected={(e, { suggestionValue }) => { this.handleAddChip(suggestionValue); e.preventDefault() }}
        focusInputOnSuggestionClick={true}
        inputProps={{
          classes,
          chips: this.getChips(),
          onChange: this.handleTextFieldInputChange,
          value: this.state.textFieldInput,
          onAdd: (chip) => this.handleAddChip(chip),
          onBeforeAdd: (chip) => this.handleBeforeAddChip(chip),
          onDelete: (chip, index) => this.handleDeleteChip(chip, index),
          ...rest
        }}
      />
    )
  }
}

export default withStyles(SearchBar.styles)(SearchBar)
