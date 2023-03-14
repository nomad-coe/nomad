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
import React, { useCallback, useState, useMemo, useRef } from 'react'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import PropTypes from 'prop-types'
import { useRecoilValue } from 'recoil'
import clsx from 'clsx'
import {
  CircularProgress,
  Tooltip,
  IconButton,
  TextField,
  ListItemText
} from '@material-ui/core'
import Autocomplete from '@material-ui/lab/Autocomplete'
import CloseIcon from '@material-ui/icons/Close'
import HelpOutlineIcon from '@material-ui/icons/HelpOutline'
import { isNil } from 'lodash'
import { useSearchContext } from '../SearchContext'
import { guiState } from '../../GUIMenu'
import { useSuggestions } from '../../../hooks'
import searchQuantities from '../../../searchQuantities'
import Placeholder from '../../visualization/Placeholder'
import { getSuggestions } from '../../../utils'

/*
 * Lowe-level representational component for all text fields used in the search.
 */
const useInputTextFieldStyles = makeStyles(theme => ({
  root: {
    height: '3rem'
  }
}))
export const InputTextField = React.memo((props) => {
  const initialLabel = useState(props.label)[0]
  const inputVariant = useRecoilValue(guiState('inputVariant'))
  const inputSize = useRecoilValue(guiState('inputSize'))
  const styles = useInputTextFieldStyles({classes: props.classes})

  return props.loading
    ? <Placeholder className={clsx(props.className, styles.root)} />
    : <TextField {...props} size={inputSize} variant={inputVariant} hiddenLabel={!initialLabel}/>
})

InputTextField.propTypes = {
  label: PropTypes.string,
  loading: PropTypes.bool,
  className: PropTypes.string,
  classes: PropTypes.object
}

/*
 * Text field that can be used to submit filter values that target a specific
 * quantity. Can also suggest values.
 */
const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    display: 'flex',
    alignItems: 'flex-start',
    justifyContent: 'center',
    flexDirection: 'column',
    boxSizing: 'border-box'
  }
}))
export const InputTextQuantity = React.memo(({
  quantity,
  suggestions,
  loading,
  onChange,
  disableSuggestions,
  className,
  classes,
  ...TextFieldProps
}) => {
  const theme = useTheme()
  const { filterData, useSetFilter } = useSearchContext()
  const styles = useStyles({classes: classes, theme: theme})
  const [inputValue, setInputValue] = useState('')
  const [suggestionInput, setSuggestionInput] = useState('')
  const [highlighted, setHighlighted] = useState({value: ''})
  const [open, setOpen] = useState(false)
  const [error, setError] = useState(false)
  const [quantitiesSuggestion, quantitiesAll] = useMemo(() => [
    [{name: quantity, size: 5}],
    new Set([quantity])
  ], [quantity])
  const [suggestionsAuto, loadingAuto] = useSuggestions(quantitiesSuggestion, quantitiesAll, suggestionInput)
  const finalSuggestions = suggestions || suggestionsAuto
  const finalLoading = loading || loadingAuto
  const disableSuggestionsFinal = suggestions
    ? true
    : isNil(disableSuggestions)
      ? !searchQuantities[quantity]?.suggestion
      : disableSuggestions

  // Attach the filter hook
  const setFilter = useSetFilter(quantity)
  const disabled = TextFieldProps.disabled

  // Sets the input value and calls the callback if given
  const handleChange = useCallback((input) => {
    setInputValue(input)
    onChange && onChange(input)
  }, [onChange])

  // Clears the input and suggestions
  const clearInputValue = useCallback(() => {
    handleChange('')
    setSuggestionInput('')
    setOpen(false)
  }, [handleChange])

  // Triggered when a value is submitted by pressing enter or clicking the
  // search icon.
  const handleSubmit = useCallback((value) => {
    if (value.trim().length === 0) {
      return
    }
    const valid = true
    if (valid) {
      // Submit to search context on successful validation.
      setFilter(old => {
        const multiple = filterData[quantity].multiple
        return (isNil(old) || !multiple) ? new Set([value]) : new Set([...old, value])
      })
      clearInputValue()
    } else {
      setError(`Invalid query`)
    }
  }, [filterData, quantity, setFilter, clearInputValue])

  const handleHighlight = useCallback((event, value, reason) => {
    setHighlighted(value)
  }, [])

  // Handles special key presses
  const handleKeyDown = useCallback((event) => {
    // When escape is pressed, close the menu if it is visible and showing some
    // items. Otherwise clear the current text input.
    if (event.key === 'Escape') {
      if (open && suggestions?.length > 0) {
        setOpen(false)
      } else {
        clearInputValue()
      }
      event.stopPropagation()
      event.preventDefault()
    }
    // When enter is pressed, select currently highlighted value and close menu,
    // or if menu is not open submit the value.
    if (event.key === 'Enter') {
      if (open && highlighted?.value) {
        handleSubmit(highlighted.value)
      } else {
        handleSubmit(inputValue)
      }
      event.stopPropagation()
      event.preventDefault()
    }
  }, [open, suggestions, highlighted, handleSubmit, inputValue, clearInputValue])

  // Handle typing events.
  const handleInputChange = useCallback((event, value, reason) => {
    setError(error => error ? undefined : null)

    // When an option is selected (mouse or keyboard), the filter is immediately
    // submitted and the field value cleared.
    if (reason === 'reset') {
      handleSubmit(value)
    } else {
      handleChange(value)
    }

    // Suggestions are only retrieved on user input, or when the value has been
    // cleared (this clears the suggestions)
    if (!disableSuggestionsFinal && (value.trim() === '' || reason === 'input')) {
      setSuggestionInput(value)
    }
  }, [disableSuggestionsFinal, handleSubmit, handleChange])

  return <div className={clsx(className, styles.root)}>
    <Autocomplete
      freeSolo
      disabled={disabled}
      clearOnBlur={false}
      inputValue={inputValue}
      value={null}
      open={open}
      onOpen={() => setOpen(true)}
      onClose={() => setOpen(false)}
      fullWidth
      disableClearable
      classes={{endAdornment: styles.endAdornment}}
      options={finalSuggestions}
      onInputChange={handleInputChange}
      onHighlightChange={handleHighlight}
      getOptionLabel={option => option.value}
      getOptionSelected={(option, value) => false}
      renderInput={(params) => {
        // We need to strip out the styling of the input field that is imposed
        // by Autocomplete. Otherwise the styles enabled by the
        // hiddenLabel-property will be overridden.
        params.InputProps.className = undefined
        return <InputTextField
          {...params}
          placeholder="Type here"
          label={error || undefined}
          error={!!error}
          onKeyDown={handleKeyDown}
          InputLabelProps={{ shrink: true }}
          InputProps={{
            ...params.InputProps,
            endAdornment: (<>
              {finalLoading ? <CircularProgress color="inherit" size={20} /> : null}
              {(inputValue?.length || null) && <>
                <Tooltip title="Clear">
                  <IconButton
                    size="small"
                    onClick={clearInputValue}
                    className={styles.iconButton}
                    aria-label="clear"
                  >
                    <CloseIcon/>
                  </IconButton>
                </Tooltip>
              </>}
            </>)
          }}
          {...TextFieldProps}
        />
      }}
    />
  </div>
})

InputTextQuantity.propTypes = {
  /*
   * The quantity targeted by the text field target.
   */
  quantity: PropTypes.string,
  /*
   * A manual array of suggestions.
   */
  suggestions: PropTypes.array,
  /*
   * Whether suggestions are being loaded.
   */
  loading: PropTypes.bool,
  /*
   * Callback for when the input changes.
   */
  onChange: PropTypes.bool,
  /*
   * Whether to enable or disable automatic suggestions. Will be forcefully set
   * to true if manual list of suggestions are provided. If no value is given
   * the suggestions are turned on if they are available for the quantity.
   */
  disableSuggestions: PropTypes.bool,
  className: PropTypes.string,
  classes: PropTypes.object
}

/*
 * Generic text field component that should be used for most user inputs.
 * Defines default behaviour for user input such as clearing inputs when
 * pressing esc and submitting values when pressing enter. Can also
 * display customizable list of suggestions.
 */
const useInputTextStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    display: 'flex',
    alignItems: 'flex-start',
    justifyContent: 'center',
    flexDirection: 'column',
    boxSizing: 'border-box'
  }
}))
export const InputText = React.memo(({
  value,
  error,
  shrink,
  suggestions,
  loading,
  onChange,
  onAccept,
  onSelect,
  onBlur,
  onError,
  disableSuggestions,
  getOptionLabel,
  renderOption,
  filterOptions,
  className,
  classes,
  ...TextFieldProps
}) => {
  const theme = useTheme()
  const styles = useInputTextStyles({classes: classes, theme: theme})
  const [open, setOpen] = useState(false)
  const disabled = TextFieldProps.disabled
  // The highlighted item is stored in a ref to keep the component more
  // responsive during browsing the suggestions
  const highlightRef = useRef(null)

  // Clears the input value and closes suggestions list
  const clearInputValue = useCallback(() => {
    onError && onError(undefined)
    onChange && onChange("")
    setOpen(false)
  }, [onChange, onError])

  // Handle item highlighting: items can he highlighted with mouse or keyboard.
  const handleHighlight = useCallback((event, value, reason) => {
    highlightRef.current = value
  }, [highlightRef])

  // Handle blur
  const handleBlur = useCallback(() => {
    onBlur && onBlur()
    onAccept && onAccept(value)
  }, [onBlur, onAccept, value])

  // Handles special key presses
  const handleKeyDown = useCallback((event) => {
    // When escape is pressed, close the menu if it is visible and showing some
    // items. Otherwise clear the current text input.
    if (event.key === 'Escape') {
      if (open && suggestions?.length > 0) {
        setOpen(false)
      } else {
        clearInputValue()
      }
      event.stopPropagation()
      event.preventDefault()
    }
    // When enter is pressed, select currently highlighted value and close menu,
    // or if menu is not open submit the value.
    if (event.key === 'Enter') {
      if (open && highlightRef.current) {
        onSelect && onSelect(getOptionLabel(highlightRef.current).trim())
      } else {
        onAccept && onAccept(value && value.trim())
      }
      event.stopPropagation()
      event.preventDefault()
      setOpen(false)
    }
  }, [open, suggestions, onSelect, onAccept, value, getOptionLabel, clearInputValue, highlightRef])

  // Handle typing events.
  const handleInputChange = useCallback((event, value, reason) => {
    onError && onError(undefined)
    event && onChange && onChange(value)
  }, [onChange, onError])

  return <div className={clsx(className, styles.root)}>
    <Autocomplete
      freeSolo
      disabled={disabled}
      clearOnBlur={false}
      inputValue={value || ''}
      value={null}
      open={open}
      onOpen={() => setOpen(true)}
      onClose={() => setOpen(false)}
      onBlur={handleBlur}
      fullWidth
      disableClearable
      classes={{endAdornment: styles.endAdornment}}
      options={suggestions}
      onInputChange={handleInputChange}
      onHighlightChange={handleHighlight}
      getOptionLabel={getOptionLabel}
      getOptionSelected={(option, value) => false}
      filterOptions={filterOptions}
      renderOption={renderOption}
      renderInput={(params) => {
        // We need to strip out the styling of the input field that is imposed
        // by Autocomplete. Otherwise the styles enabled by the
        // hiddenLabel-property will be overridden.
        params.InputProps.className = undefined
        return <InputTextField
          {...params}
          helperText={error || undefined}
          error={!!error}
          onKeyDown={handleKeyDown}
          InputLabelProps={{shrink}}
          InputProps={{
            ...params.InputProps,
            endAdornment: (<>
              {loading ? <CircularProgress color="inherit" size={20} /> : null}
              {(value?.length || null) && <>
                <Tooltip title="Clear">
                  <IconButton
                    size="small"
                    onClick={clearInputValue}
                    className={styles.iconButton}
                    aria-label="clear"
                  >
                    <CloseIcon/>
                  </IconButton>
                </Tooltip>
              </>}
            </>)
          }}
          {...TextFieldProps}
        />
      }}
    />
  </div>
})

InputText.propTypes = {
  value: PropTypes.string,
  error: PropTypes.string, // Error shown underneath the text
  shrink: PropTypes.bool, // Whether the label should automatically "shrink" on input
  suggestions: PropTypes.array, // Array of suggested values
  loading: PropTypes.bool, // Whether loading icon should be shown
  onChange: PropTypes.func, // Triggered whenever the input text changes
  onSelect: PropTypes.func, // Triggered when an option is selected from suggestions
  onAccept: PropTypes.func, // Triggered when value should be accepted
  onBlur: PropTypes.func, // Triggered when text goes out of focus
  onError: PropTypes.func, // Triggered when any errors should be cleared
  disableSuggestions: PropTypes.bool,
  getOptionLabel: PropTypes.func,
  renderOption: PropTypes.func,
  filterOptions: PropTypes.func,
  className: PropTypes.string,
  classes: PropTypes.object
}

InputText.defaultProps = {
  getOptionLabel: (option) => option.value
}

/**
 * Wrapper around InputMetainfo which automatically shows suggestions and only
 * accepts metainfo that exist in the current search context.
 */
export const InputSearchMetainfo = React.memo(({
  value,
  label,
  error,
  onChange,
  onSelect,
  onAccept,
  onError,
  dtypes,
  dtypesRepeatable,
  optional,
  disableNonAggregatable
}) => {
  const { filterData } = useSearchContext()

  // Fetch the available metainfo names and create options that are compatible
  // with InputMetainfo.
  const suggestions = useMemo(() => {
    const suggestions = Object.entries(filterData)
      .filter(([key, data]) => {
        if (disableNonAggregatable && !data.aggregatable) return false
        const dtype = data?.dtype
        return data?.repeatsRecursive
          ? dtypesRepeatable?.has(dtype)
          : dtypes?.has(dtype)
      })
      .map(([key, data]) => ({path: key, description: data.description}))
    return suggestions
  }, [filterData, dtypes, dtypesRepeatable, disableNonAggregatable])

  return <InputMetainfo
    options={suggestions}
    value={value}
    label={label}
    error={error}
    onChange={onChange}
    onSelect={onSelect}
    onAccept={onAccept}
    onError={onError}
    optional={optional}
  />
})

InputSearchMetainfo.propTypes = {
  label: PropTypes.string,
  value: PropTypes.string,
  error: PropTypes.string,
  onChange: PropTypes.func,
  onSelect: PropTypes.func,
  onAccept: PropTypes.func,
  onError: PropTypes.func,
  /* List of allowed data types for non-repeatable quantities. */
  dtypes: PropTypes.object,
  /* List of allowed data types for repeatable quantities. */
  dtypesRepeatable: PropTypes.object,
  /* Whether the value is optional */
  optional: PropTypes.bool,
  /* Whether non-aggregatable values are excluded */
  disableNonAggregatable: PropTypes.bool
}

/**
 * Wrapper around InputText that is specialized in showing metainfo options.
 */
const itemMargin = 6
export const useInputStyles = makeStyles(theme => ({
  optionText: {
    flexGrow: 1
  },
  path: {
    // Allows very long metainfo names to break into several lines
    wordBreak: 'break-all'
  },
  option: {
    // The negative margins ensure that the item covers the entire container
    // size, which ensures that hover is always active when mouse is inside it.
    marginTop: -itemMargin,
    marginBottom: -itemMargin,
    paddingTop: itemMargin,
    paddingBottom: itemMargin,
    width: '100%',
    display: 'flex',
    alignItems: 'stretch',
    // The description icon is hidden until the item is hovered. It is not
    // removed from the document with "display: none" in order for the hover to
    // not change the layout which may cause other elements to shift around.
    '& .description': {
      visibility: "hidden",
      display: 'flex',
      width: theme.spacing(5),
      marginLeft: theme.spacing(1),
      marginTop: -itemMargin,
      marginBottom: -itemMargin,
      alignItems: 'center',
      justifyContent: 'center'
    },
    '&:hover .description': {
      visibility: "visible"
    }
  }
}))
export const InputMetainfo = React.memo(({
  label,
  value,
  options,
  error,
  onChange,
  onSelect,
  onAccept,
  onBlur,
  onError,
  optional
}) => {
  const styles = useInputStyles()

  // Predefine all option objects, all option paths and also pre-tokenize the
  // options for faster matching.
  const { optionsMap, paths, pathsSet, filter } = useMemo(() => {
    const optionsMap = Object.fromEntries(options.map((option) => {
      return [option.path, {...option, path: option.path}]
    }))
    const paths = Object.keys(optionsMap)
    const pathsSet = new Set(paths)
    const { filter } = getSuggestions(paths, 0)
    return { optionsMap, paths, pathsSet, filter }
  }, [options])

  // Used to validate the input and raise errors
  const validate = useCallback((value) => {
    const empty = !value || value.length === 0
    if (optional && empty) {
      return {valid: true, error: undefined}
    } else if (empty) {
      return {valid: false, error: 'Please specify a value.'}
    } else if (!(pathsSet.has(value))) {
      return {valid: false, error: 'Invalid value for this field.'}
    }
    return {valid: true, error: undefined}
  }, [pathsSet, optional])

  // Handles the final acceptance of a value
  const handleAccept = useCallback((value) => {
    const {valid, error} = validate(value)
    if (valid) {
      onAccept && onAccept(value)
    } else {
      onError && onError(error)
    }
  }, [validate, onError, onAccept])

  return <InputText
    value={value || null}
    label={label}
    error={error}
    onChange={onChange}
    onSelect={onSelect}
    onAccept={handleAccept}
    onBlur={onBlur}
    onError={onError}
    suggestions={paths}
    getOptionLabel={option => option}
    filterOptions={(options, { inputValue }) => filter(inputValue).map(option => option.value)}
    renderOption={(path) => {
      const option = optionsMap[path]
      return <div className={styles.option}>
        <ListItemText
          primary={option.path}
          secondary={option.secondary}
          className={styles.optionText}
          primaryTypographyProps={{className: styles.path}}
        />
        {option.description &&
          <Tooltip title={option.description || ''}>
            <div className="description">
              <HelpOutlineIcon fontSize="small" color="action"/>
            </div>
          </Tooltip>
        }
      </div>
    }}
  />
})

InputMetainfo.propTypes = {
  label: PropTypes.string,
  value: PropTypes.string,
  options: PropTypes.arrayOf(PropTypes.shape({
    path: PropTypes.string,
    secondary: PropTypes.string,
    description: PropTypes.string
  })),
  error: PropTypes.string,
  onChange: PropTypes.func,
  onSelect: PropTypes.func,
  onAccept: PropTypes.func,
  onError: PropTypes.func,
  onBlur: PropTypes.func,
  optional: PropTypes.bool // Set to true if field can be empty
}

InputMetainfo.defaultProps = {
  label: "quantity"
}
