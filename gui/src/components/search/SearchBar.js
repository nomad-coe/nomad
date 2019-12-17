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
import { nomadPrimaryColor } from '../../config'
import { withDomain } from '../domains'
import { compose } from 'recompose'
import SearchContext from '../search/SearchContext'
import { withApi } from '../api'


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
              backgroundColor: isFocused ? nomadPrimaryColor[500] : undefined,
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
    domain: PropTypes.object.isRequired,
    loading: PropTypes.number
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

  constructor(props) {
    super(props)
    const { domain } = props
    const reStr = `^(${Object.keys(domain.additionalSearchKeys).join('|')})=`
    this.additionalSearchKeyRE = new RegExp(reStr)
  }

  getSuggestions(valueWithCase) {
    const value = valueWithCase.toLowerCase()

    const {statistics} = this.context.state.response
    const suggestions = []

    // filter out pseudo quantity total
    const quantityKeys = Object.keys(statistics).filter(quantity => quantity !== 'total')

    // put authors to the end
    const authorIndex = quantityKeys.indexOf('authors')
    if (authorIndex >= 0) {
      quantityKeys[authorIndex] = quantityKeys.splice(quantityKeys.length - 1, 1, quantityKeys[authorIndex])[0]
    }

    quantityKeys.forEach(quantity => {
      Object.keys(statistics[quantity]).forEach(quantityValue => {
        const quantityValueLower = quantityValue.toLowerCase()
        if (quantityValueLower.startsWith(value) || (quantity === 'authors' && quantityValueLower.includes(value))) {
          suggestions.push({
            key: quantity,
            value: quantityValue
          })
        }
      })
    })

    // Add additional quantities to the end
    const match = value.match(this.additionalSearchKeyRE)
    if (match && this.props.domain.additionalSearchKeys[match[1]]) {
      suggestions.push({
        key: match[1],
        value: valueWithCase.substring(match[0].length)
      })
    }

    // Always add as comment to the end of suggestions
    suggestions.push({
      key: 'comment',
      value: value
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
    const values = {...this.context.state.query}

    let key, value
    if (chip.includes('=')) {
      const parts = chip.split(/=(.+)/)
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
    this.context.setQuery(values, true)
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

    const {state: {query}, setQuery} = this.context
    const values = {...query}
    delete values[key]
    setQuery(values, true)
  }

  getChips() {
    const {state: {query: {owner, ...values}}} = this.context
    return Object.keys(values).filter(key => values[key]).map(key => {
      if (key === 'atoms') {
        return `atoms=[${values[key].join(',')}]`
      } else {
        return `${key}=${values[key]}`
      }
    })
  }

  static contextType = SearchContext.type

  render() {
    const {classes, domain, loading} = this.props
    const {response: {pagination, statistics}, query} = this.context.state

    let helperText = <span>loading ...</span>
    if (pagination && statistics) {
      if (pagination.total === 0) {
        helperText = <span>There are no more entries matching your criteria.</span>
      } else {
        helperText = <span>
          There {pagination.total === 1 ? 'is' : 'are'} {Object.keys(domain.searchMetrics).filter(key => statistics.total.all[key]).map(key => {
            return <span key={key}>
              {domain.searchMetrics[key].renderResultString(!loading ? statistics.total.all[key] : '...')}
            </span>
          })}{Object.keys(query).length ? ' left' : ''}.
        </span>
      }
    }

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
          label: 'search',
          fullWidth: true,
          fullWidthInput: false,
          InputLabelProps: {
            shrink: true
          },
          placeholder: domain.searchPlaceholder,
          helperText: helperText
        }}
      />
    )
  }
}

export default compose(withApi(false, false), withDomain, withStyles(SearchBar.styles))(SearchBar)
