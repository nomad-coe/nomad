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
import React, { useCallback, useState, useMemo, useRef, useEffect } from 'react'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { CircularProgress, Tooltip, IconButton, TextField } from '@material-ui/core'
import Autocomplete from '@material-ui/lab/Autocomplete'
import ArrowDropDownIcon from '@material-ui/icons/ArrowDropDown'
import CloseIcon from '@material-ui/icons/Close'
import { isNil } from 'lodash'
import { useSearchContext } from '../SearchContext'
import { useSuggestions } from '../../../hooks'
import { searchQuantities } from '../../../config'
import Placeholder from '../../visualization/Placeholder'

/*
 * Low-level representational component for all text fields used in the search.
 */
const useInputTextFieldStyles = makeStyles(theme => ({
  root: {
    height: '3rem'
  }
}))
export const InputTextField = React.memo((props) => {
  const initialLabel = useState(props.label)[0]
  const styles = useInputTextFieldStyles({classes: props.classes})

  return props.loading
    ? <Placeholder className={clsx(props.className, styles.root)} />
    : <TextField size="small" variant="filled" {...props} hiddenLabel={!initialLabel}/>
})

InputTextField.propTypes = {
  label: PropTypes.string,
  loading: PropTypes.bool,
  className: PropTypes.string,
  classes: PropTypes.object
}

/*
 * Customized version of Autocomplete with custom NOMAD styling and behaviour.
 *
 * Defines default behaviour for user input such as clearing inputs when
 * pressing esc and submitting values when pressing enter. Can also display
 * customizable list of suggestions.
 */
const useInputTextStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    display: 'flex',
    alignItems: 'flex-start',
    justifyContent: 'center',
    flexDirection: 'column',
    boxSizing: 'border-box'
  },
  popupIndicatorOpen: {
    transform: 'rotate(180deg)'
  },
  adornmentList: {
    display: 'flex',
    alignItems: 'center'
  },
  adornment: {
    padding: '3px'
  },
  listbox: {
    boxSizing: 'border-box',
    '& ul': {
      padding: 0,
      margin: 0
    }
  },
  option: {
    paddingTop: 0,
    paddingBottom: 0
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
  onHighlight,
  onBlur,
  onFocus,
  onError,
  getOptionLabel,
  groupBy,
  renderOption,
  renderGroup,
  suggestAllOnFocus,
  showOpenSuggestions,
  autoHighlight,
  ListboxComponent,
  filterOptions,
  className,
  classes,
  TextFieldProps,
  InputProps,
  PaperComponent,
  disableClearable,
  disableAcceptOnBlur,
  validate,
  disableValidateOnSelect
}) => {
  const theme = useTheme()
  const styles = useInputTextStyles({classes: classes, theme: theme})
  const [open, setOpen] = useState(false)
  const [suggestAll, setSuggestAll] = useState(false)
  const disabled = TextFieldProps?.disabled
  // The highlighted item is stored in a ref to keep the component more
  // responsive during browsing the suggestions
  const highlightRef = useRef(null)

  // Clears the input value and closes suggestions list
  const clearInputValue = useCallback(() => {
    onError?.(undefined)
    onChange?.("")
    setOpen(false)
  }, [onChange, onError])

  const handleAccept = useCallback((value) => {
    const {valid, error, data} = validate ? validate(value) : {valid: true, error: undefined, data: undefined}
    if (valid) {
      onAccept?.(value && value.trim(), data)
    } else {
      onError?.(error)
    }
  }, [onAccept, validate, onError])

  // Validate the initial value if it is non-empty.
  useEffect(() => {
    if (!(isNil(value) || value?.trim?.() === '')) {
      handleAccept(value)
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  const handleSelect = useCallback((value) => {
    if (disableValidateOnSelect) {
      onSelect?.(value && value.trim())
    } else {
      const {valid, error, data} = validate ? validate(value) : {valid: true, error: undefined, data: undefined}
      if (valid) {
        onSelect?.(value && value.trim(), data)
      } else {
        onError?.(error)
      }
    }
  }, [onSelect, disableValidateOnSelect, validate, onError])

  // Handle item highlighting: items can he highlighted with mouse or keyboard.
  const handleHighlight = useCallback((event, value, reason) => {
    onHighlight?.(value, reason)
    highlightRef.current = value
  }, [highlightRef, onHighlight])

  // Handle blur
  const handleBlur = useCallback(() => {
    onBlur?.()
    !disableAcceptOnBlur && handleAccept(value)
  }, [onBlur, handleAccept, value, disableAcceptOnBlur])

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
        handleSelect?.(getOptionLabel(highlightRef.current).trim())
      } else {
        handleAccept(value && value.trim())
      }
      event.stopPropagation()
      event.preventDefault()
      setOpen(false)
    }
  }, [open, suggestions, handleSelect, handleAccept, value, getOptionLabel, clearInputValue, highlightRef])

  // Handle input events. Errors are cleaned in input change, regular typing
  // emits onChange, selection with mouse emits onSelect.
  const handleInputChange = useCallback((event, value, reason) => {
    setSuggestAll(false)
    onError && onError(undefined)
    if (event) {
      if (reason === 'reset') {
        handleSelect?.(value)
      } else {
        onChange?.(value)
      }
    }
  }, [onChange, handleSelect, onError])

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
      autoHighlight={autoHighlight}
      classes={{
        listbox: styles.listbox,
        option: styles.option
      }}
      ListboxComponent={ListboxComponent}
      PaperComponent={PaperComponent}
      options={suggestions}
      onInputChange={handleInputChange}
      onHighlightChange={handleHighlight}
      getOptionLabel={getOptionLabel}
      getOptionSelected={(option, value) => false}
      groupBy={groupBy}
      renderGroup={renderGroup}
      filterOptions={suggestAll
          ? (opt) => opt
          : filterOptions
      }
      renderOption={renderOption}
      selectOnFocus={true}
      renderInput={(params) => {
        // We need to strip out the styling of the input field that is imposed
        // by Autocomplete. Otherwise the styles enabled by the
        // hiddenLabel-property will be overridden.
        params.InputProps.className = undefined
        return <InputTextField
          {...params}
          size="small"
          helperText={error || undefined}
          error={!!error}
          onFocus={() => { suggestAllOnFocus && setSuggestAll(true); onFocus?.() } }
          onKeyDown={handleKeyDown}
          InputLabelProps={{shrink}}
          InputProps={{
            ...params.InputProps,
            ...InputProps,
            endAdornment: (<div className={styles.adornmentList}>
              {loading ? <CircularProgress color="inherit" size={20} className={styles.adornment} /> : null}
              {(value?.length && !disableClearable)
                ? <Tooltip title="Clear">
                    <IconButton
                      size="small"
                      disabled={disabled}
                      onClick={clearInputValue}
                      className={styles.iconButton}
                      aria-label="clear"
                    >
                      <CloseIcon/>
                    </IconButton>
                  </Tooltip>
                : null
              }
              {(showOpenSuggestions)
                ? <IconButton
                    size="small"
                    disabled={disabled}
                    onClick={() => setOpen(old => !old)}
                    className={clsx(styles.popupIndicator, {
                      [styles.popupIndicatorOpen]: open
                    })}
                  >
                    <ArrowDropDownIcon />
                </IconButton>
                : null
              }
              {InputProps?.endAdornment || null}
            </div>)
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
  onFocus: PropTypes.func, // Triggered when text is focused
  onHighlight: PropTypes.func, // Triggered when selection is highlighted
  onError: PropTypes.func, // Triggered when any errors should be cleared
  getOptionLabel: PropTypes.func,
  groupBy: PropTypes.func,
  renderOption: PropTypes.func,
  renderGroup: PropTypes.func,
  ListboxComponent: PropTypes.any,
  PaperComponent: PropTypes.any,
  TextFieldProps: PropTypes.object,
  InputProps: PropTypes.object,
  filterOptions: PropTypes.func,
  disableClearable: PropTypes.bool,
  autoHighlight: PropTypes.bool,
  disableAcceptOnBlur: PropTypes.bool,
  validate: PropTypes.func, // Function that can be used to validate the input
  disableValidateOnSelect: PropTypes.bool, // Whether validation on selecting autocompletion value should be disabled
  suggestAllOnFocus: PropTypes.bool, // Whether to provide all suggestion values when input is focused
  showOpenSuggestions: PropTypes.bool, // Whether to show button for opening suggestions
  className: PropTypes.string,
  classes: PropTypes.object
}

InputText.defaultProps = {
  getOptionLabel: (option) => option.value,
  showOpenSuggestions: false,
  autoHighlight: false
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
  const [suggestionsAuto, loadingAuto] = useSuggestions(quantitiesSuggestion, quantitiesAll, suggestionInput, filterData)
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
          size="small"
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
