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
import React, { useCallback, useState, useMemo } from 'react'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import PropTypes from 'prop-types'
import { useRecoilValue } from 'recoil'
import clsx from 'clsx'
import {
  CircularProgress,
  Tooltip,
  IconButton,
  TextField
} from '@material-ui/core'
import Autocomplete from '@material-ui/lab/Autocomplete'
import CloseIcon from '@material-ui/icons/Close'
import { isNil } from 'lodash'
import { useSearchContext } from '../SearchContext'
import { guiState } from '../../GUIMenu'
import { useSuggestions } from '../../../hooks'
import searchQuantities from '../../../searchQuantities'

/*
 * Representational component for all text fields used in the search.
 */
export const InputTextField = React.memo((props) => {
  const initialLabel = useState(props.label)[0]
  const inputVariant = useRecoilValue(guiState('inputVariant'))
  const inputSize = useRecoilValue(guiState('inputSize'))
  return <TextField {...props} size={inputSize} variant={inputVariant} hiddenLabel={!initialLabel}/>
})

InputTextField.propTypes = {
  label: PropTypes.string
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
  const { filterData, useSetFilter, useFilterLocked } = useSearchContext()
  const styles = useStyles({classes: classes, theme: theme})
  const [inputValue, setInputValue] = useState('')
  const [suggestionInput, setSuggestionInput] = useState('')
  const [highlighted, setHighlighted] = useState({value: ''})
  const [open, setOpen] = useState(false)
  const [error, setError] = useState(false)
  const suggestionQuantity = useMemo(() => [quantity], [quantity])
  const [suggestionsAuto, loadingAuto] = useSuggestions(suggestionQuantity, suggestionInput)
  const finalSuggestions = suggestions || suggestionsAuto
  const finalLoading = loading || loadingAuto
  const disableSuggestionsFinal = suggestions
    ? true
    : isNil(disableSuggestions)
      ? !searchQuantities[quantity]?.suggestion
      : disableSuggestions

  // Attach the filter hook
  const setFilter = useSetFilter(quantity)
  const locked = useFilterLocked(quantity)
  const disabled = locked || TextFieldProps.disabled

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
    let valid = true
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
   * Whether suggestions are being loade
   */
  loading: PropTypes.bool,
  /*
   * Callback for when the input changes
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
