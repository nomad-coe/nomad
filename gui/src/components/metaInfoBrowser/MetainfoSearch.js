import React from 'react'
import PropTypes from 'prop-types'
import Autosuggest from 'react-autosuggest'
import parse from 'autosuggest-highlight/parse'
import TextField from '@material-ui/core/TextField'
import Paper from '@material-ui/core/Paper'
import MenuItem from '@material-ui/core/MenuItem'
import { withStyles } from '@material-ui/core/styles'

function renderInputComponent(inputProps) {
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

function renderSuggestion(suggestion, { query, isHighlighted }) {
  const inputValue = query.trim()
  const matches = match(suggestion.name, inputValue)
  const parts = parse(suggestion.name, matches)

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

function getSuggestionValue(suggestion) {
  return suggestion.name
}

const styles = theme => ({
  container: {
    position: 'relative'
  },
  suggestionsContainerOpen: {
    position: 'absolute',
    zIndex: 1000,
    marginTop: theme.spacing(1),
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
  divider: {
    height: theme.spacing(2)
  },
  input: {
    height: 24
  }
})

function match(content, query) {
  const queries = query.toLowerCase().split(' ')
  const result = queries.map(query => {
    const index = content.toLowerCase().indexOf(query)
    if (index >= 0) {
      return [index, index + query.length]
    } else {
      return null
    }
  }).filter(match => match !== null)

  if (result.length === queries.length) {
    result.sort((a, b) => a[0] - b[0])
    return result
  } else {
    return []
  }
}

class MetainfoSearch extends React.Component {
  state = {
    single: '',
    popper: '',
    suggestions: []
  }

  lastRequestTimeout = null

  getSuggestions(value) {
    const inputValue = value.trim()
    return this.props.suggestions.filter(suggestion =>
      match(suggestion.name, inputValue).length > 0
    ).slice(0, 15)
  }

  handleSuggestionsFetchRequested = ({ value }) => {
    if (this.lastRequestTimeout !== null) {
      clearTimeout(this.lastRequestTimeout)
    }

    this.lastRequestTimeout = setTimeout(() => {
      this.setState({
        suggestions: this.getSuggestions(value)
      })
    }, 200)
  }

  handleSuggestionsClearRequested = () => {
    this.setState({
      suggestions: []
    })
  }

  handleChange = name => (event, { newValue }) => {
    this.setState({
      [name]: newValue
    })
    this.props.onChange(newValue)
  }

  render() {
    const { classes } = this.props

    const autosuggestProps = {
      renderInputComponent,
      suggestions: this.state.suggestions,
      onSuggestionsFetchRequested: this.handleSuggestionsFetchRequested,
      onSuggestionsClearRequested: this.handleSuggestionsClearRequested,
      getSuggestionValue,
      renderSuggestion
    }

    return (
      <Autosuggest
        className={classes.root}
        {...autosuggestProps}
        inputProps={{
          classes: classes,
          label: 'Definition',
          placeholder: 'Enter definition name',
          value: this.state.single,
          onChange: this.handleChange('single'),
          InputLabelProps: {
            shrink: true
          },
          fullWidth: true
        }}
        theme={{
          container: classes.container,
          suggestionsContainerOpen: classes.suggestionsContainerOpen,
          suggestionsList: classes.suggestionsList,
          suggestion: classes.suggestion
        }}
        renderSuggestionsContainer={options => (
          <Paper {...options.containerProps} square>
            {options.children}
          </Paper>
        )}
      />
    )
  }
}

MetainfoSearch.propTypes = {
  classes: PropTypes.object.isRequired,
  suggestions: PropTypes.array.isRequired,
  onChange: PropTypes.func.isRequired
}

export default withStyles(styles)(MetainfoSearch)
