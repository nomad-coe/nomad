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
import React, { useCallback, useState, useMemo, createContext, useContext, useRef } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import assert from 'assert'
import { isNil, has } from 'lodash'
import { makeStyles } from '@material-ui/core/styles'
import SearchIcon from '@material-ui/icons/Search'
import HistoryIcon from '@material-ui/icons/History'
import CloseIcon from '@material-ui/icons/Close'
import { HelpButton } from '../../components/Help'
import { Paper, Tooltip, Chip, List, ListSubheader, Typography, IconButton } from '@material-ui/core'
import { parseQuantityName, getSchemaAbbreviation } from '../../utils'
import { useSuggestions } from '../../hooks'
import { useSearchContext } from './SearchContext'
import { quantityNameSearch } from './FilterRegistry'
import { VariableSizeList } from 'react-window'
import { InputText } from './input/InputText'
import { SearchSuggestion, SuggestionType } from './SearchSuggestion'
import { SearchSyntaxes } from './SearchSyntaxes'
import { renderRow, useResetCache, LISTBOX_PADDING } from './input/InputMetainfo'

/**
 * Displays a suggestion in the search bar.
 */
export const useSuggestionStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    height: '100%',
    display: 'flex',
    alignItems: 'center',
    '& > *': {
        marginRight: theme.spacing(1)
    },
    '& > *:last-child': {
        marginRight: 0
    },
    // The delete icon is hidden until the item is hovered. It is not removed
    // from the document with "display: none" in order for the hover to not
    // change the layout which may cause other elements to shift around.
    '& .delete': {
      visibility: "hidden"
    },
    '&:hover .delete': {
      visibility: "visible"
    }
  },
  primary: {
    flexGrow: 1,
    whiteSpace: 'nowrap',
    overflowX: 'scroll',
    '&::-webkit-scrollbar': {
      display: 'none'
    },
    '-ms-overflow-style': 'none',
    scrollbarWidth: 'none'
  }
}))
const Suggestion = ({suggestion, onDelete, tooltip}) => {
  const typeTooltip = (suggestion.type === SuggestionType.Name
    ? 'Suggested name for an available filter in this app'
    : SearchSyntaxes[suggestion.type]?.readme) || ''
  const label = (suggestion.type === SuggestionType.Name
    ? 'name'
    : SearchSyntaxes[suggestion.type]?.labelShort) || ''
  const styles = useSuggestionStyles()
  let schema, primary
  if (suggestion.type === SuggestionType.Name) {
    const {path, schema: schemaTmp} = parseQuantityName(suggestion.input)
    primary = path
    schema = schemaTmp
  } else {
    primary = suggestion.input
  }

  return <div className={styles.root}>
    <Typography className={styles.primary}>
      <Tooltip title={tooltip || ''} placement='bottom' enterDelay={500} enterNextDelay={500}>
        <span>{primary}</span>
      </Tooltip>
    </Typography>
    {suggestion.history
      ? <Tooltip title="Remove from search history">
        <IconButton className="delete" size="small" onClick={(event) => {
          event.preventDefault()
          event.stopPropagation()
          onDelete?.()
        }}>
          <CloseIcon fontSize="small" color="action" />
        </IconButton>
      </Tooltip>
      : null
    }
    {suggestion.history
      ? <Tooltip title="From search history">
        <HistoryIcon fontSize="small" color="action"/>
      </Tooltip>
      : null
    }
    {schema
      ? <Tooltip title={`Definition from ${schema}`}>
      <Typography variant='caption'>{getSchemaAbbreviation(schema)}</Typography>
    </Tooltip>
      : null
    }
    <Tooltip title={typeTooltip}>
      <Chip size='small' label={label}/>
    </Tooltip>
  </div>
}

Suggestion.propTypes = {
  suggestion: PropTypes.object,
  onDelete: PropTypes.func,
  tooltip: PropTypes.string
}

export const useStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    alignItems: 'center',
    position: 'relative'
  },
  iconButton: {
    padding: 10
  }
}))

/**
 * This component shows a search bar with autocomplete functionality.
 */
const SearchBar = React.memo(({
  quantities,
  className
}) => {
  const styles = useStyles()
  const {
    filters,
    filterData,
    filterFullnames,
    useUpdateFilter,
    useParseQuery,
    useSearchSuggestions,
    usePushSearchSuggestion,
    useRemoveSearchSuggestion,
    useSetPagination,
    searchSyntaxes
  } = useSearchContext()
  const includedFormats = Object
    .keys(SearchSyntaxes)
    .filter((key) => !searchSyntaxes?.exclude?.includes(key))
  const formatReadmeList = includedFormats
    .map((key) => {
      const examples = SearchSyntaxes[key]?.examples?.map((example) => `   ${example}`).join('\n')
      return ` - **${SearchSyntaxes[key]?.label}**: ${SearchSyntaxes[key]?.readme} For example:\n\n   \`\`\`\n${examples}\n   \`\`\`\n`
    })
    .join('\n')
  const [inputValue, setInputValue] = useState('')
  const [suggestionInput, setSuggestionInput] = useState('')
  const [error, setError] = useState(false)
  const setPagination = useSetPagination()
  const updateFilter = useUpdateFilter()
  const parseQuery = useParseQuery()
  const inputRef = useRef()

  const [quantitiesAll, quantitiesAllSet, quantitiesSuggestable] = useMemo(() => {
    let quantitySetSuggestable, quantityList
    // Custom list of quantities. They are validated against the current search
    // context.
    if (quantities) {
      for (const quantity of quantities) {
        assert(
          quantity.name === quantityNameSearch || has(filterData, quantity.name),
          `Quantity ${quantity.name} does not exist in the current search context.`
        )
      }
      quantityList = [...quantities]
      quantitySetSuggestable = new Set(
        quantityList
          .map(q => q.name)
          .filter(name => name === quantityNameSearch || filterData[name]?.suggestion)
      )
    // Default quantities to use. Certain quantities are prioritized, others
    // come in alphabetical order.
    } else {
      quantitySetSuggestable = new Set([
        'results.material.elements',
        'results.material.chemical_formula_hill',
        'results.material.chemical_formula_anonymous',
        'results.material.structural_type',
        'results.material.symmetry.structure_name',
        'results.method.simulation.program_name',
        'authors.name'
      ])
      const filterList = [...filters]
      filterList.forEach(q => quantitySetSuggestable.add(q))
      quantityList = filterList.map((name) => ({name, size: 5}))
      // The list of available quantity names is provided at the very
      // bottom.
      quantitySetSuggestable.add(quantityNameSearch)
      quantityList.push({name: quantityNameSearch})
    }
    const quantitySet = new Set(quantityList.map((q) => q.name))
    return [quantityList, quantitySet, quantitySetSuggestable]
  }, [quantities, filterData, filters])

  const [suggestionNames, setSuggestionNames] = useState(quantitiesSuggestable)
  const suggestionQuantities = useMemo(() => {
    return quantitiesAll.filter((q) => suggestionNames.has(q.name))
  }, [quantitiesAll, suggestionNames])

  const [suggestionsMatch, setSuggestionsMatch] = useState([])
  const pushSuggestion = usePushSearchSuggestion()
  const removeSuggestion = useRemoveSearchSuggestion()
  const suggestionsHistory = useSearchSuggestions(inputValue)
  const [suggestionsMetainfo, loading] = useSuggestions(suggestionQuantities, quantitiesAllSet, suggestionInput, filterData)

  // Contains the final set of suggestions
  const [suggestions, keys] = useMemo(() => {
    const suggestions = {}

    for (const suggestion of suggestionsMatch) {
      suggestions[suggestion.key] = suggestion
    }

    for (const suggestion of suggestionsHistory) {
      suggestions[suggestion.key] = suggestion
    }

    for (const suggestion of suggestionsMetainfo) {
      let input, type
      if (suggestion.category === quantityNameSearch) {
        type = SuggestionType.Name
        input = suggestion.value
      } else {
        type = SuggestionType.Equality
        input = `${suggestion.category} = ${suggestion.value}`
      }
      const suggestionObject = new SearchSuggestion({input, type})
      const key = suggestionObject.key
      if (!has(suggestions, key)) {
        suggestions[key] = suggestionObject
      }
    }
    return [suggestions, Object.keys(suggestions)]
  }, [suggestionsMatch, suggestionsHistory, suggestionsMetainfo])

  const clearSuggestions = useCallback(() => {
    setSuggestionInput('')
    setSuggestionsMatch([])
  }, [])

  const clearInput = useCallback(() => {
    setError(error => error ? undefined : null)
    setInputValue('')
    setSuggestionInput('')
  }, [])

  // Used to check the validity of the given quantity name
  const checkMetainfo = useCallback((name) => {
    const fullName = filterFullnames[name] || name
    if (!quantitiesAllSet.has(fullName)) {
      return [undefined, `The quantity '${name}' is not supported`]
    }
    if (filterData[fullName].section) {
      return [fullName, `Cannot target metainfo sections`]
    }
    return [fullName, undefined]
  }, [quantitiesAllSet, filterData, filterFullnames])

  // Get suggestions based on the current input and the supported query syntaxes
  const getSuggestionsMatch = useCallback((input) => {
    input = input?.trim()
    if (!input?.length) {
      return []
    }
    const suggestions = []
    for (const key of includedFormats) {
        // regex.test() is significantly faster than string.match()
        const match = SearchSyntaxes[key].regex.test(input)
        if (match) {
          suggestions.push(new SearchSuggestion({input, type: key}))
        }
    }
    return suggestions
  }, [includedFormats])

  // Triggered when the user selects a highlighted item with mouse or with
  // keyboard.
  const handleAccept = useCallback((key) => {
    let suggestion = suggestions[key]
    if (!suggestion) {
      const bestSuggestion = suggestions[keys[0]]
      if (bestSuggestion?.input?.trim() === key?.trim()) {
        suggestion = bestSuggestion
      } else {
        setError('Unsupported query')
        return
      }
    }
    if (suggestion.type === SuggestionType.Name) {
      setInputValue(suggestion.input)
      return
    }

    // Parse the suggestion using the matched format
    const format = SearchSyntaxes[suggestion.type]
    const {target: targetTmp, value: valueTmp, inputNormalized} = format.parse(suggestion.input)

    // Check that targeted metainfo is ok, and attempt to parse the query
    const [target, errorMeta] = checkMetainfo(targetTmp)
    let value, errorQuery
    if (!errorMeta) {
      try {
        value = parseQuery(target, valueTmp)
      } catch (err) {
        errorQuery = `Invalid value for the metainfo ${target}. Please check your syntax.`
      }
    }

    // Catch any errors and raise them
    const error = errorMeta || errorQuery
    if (error) {
      setError(error)
      return
    }

    // The pagination is explicitly set to be done by descending relevance
    // score if a free-text query is performed
    if (suggestion.type === SuggestionType.Freetext) {
      setPagination(old => {
        return {...old, order_by: '_score', order: 'desc'}
      })
    }

    // Update the search filter
    updateFilter([target, old => {
      const multiple = filterData[target].multiple
      return (format?.label === SuggestionType.RangeHalfBounded || suggestion.type === SuggestionType.RangeBounded || isNil(old) || !multiple)
        ? value
        : new Set([...old, ...value])
    }])

    // We need to create a copy here, since Recoil won't allow mutating object
    // stored in a state.
    const normalizedSuggestion = new SearchSuggestion(suggestion)
    normalizedSuggestion.input = inputNormalized

    pushSuggestion(normalizedSuggestion)
    clearSuggestions()
    clearInput()
  }, [checkMetainfo, clearInput, clearSuggestions, filterData, keys, parseQuery, pushSuggestion, suggestions, updateFilter, setPagination])

  // When a suggestion is highlighted with keyboard (not mouse), the text field
  // is changed to that suggestion value. Note that this should not trigger the
  // suggestions etc. and this is why we don't set inputValue.
  const handleHighlight = useCallback((key, reason) => {
    if (reason === 'keyboard') {
      const suggestion = suggestions[key]?.input
      if (suggestion) {
        inputRef.current.value = suggestion
      }
    }
  }, [suggestions])

  // Handle typing events
  const handleInputChange = useCallback((input) => {
    setError(error => error ? undefined : null)
    setInputValue(input)
    input = input?.trim()
    if (!input) {
      clearSuggestions()
      return
    }

    // Get suggestions for making specific queries
    setSuggestionsMatch(getSuggestionsMatch(input))

    // Start getting suggestions from ES
    const split = input.split('=', 2)
    let quantitySet = quantitiesSuggestable
    if (split.length === 2) {
      const quantityName = split[0].trim()
      const quantityFullname = filterFullnames[quantityName]
      if (quantitiesSuggestable.has(quantityName)) {
        quantitySet = new Set([quantityName])
        input = split[1].trim()
      } else if (quantitiesSuggestable.has(quantityFullname)) {
        quantitySet = new Set([quantityFullname])
        input = split[1].trim()
      }
    }
    setSuggestionNames(quantitySet)
    setSuggestionInput(input)
  }, [quantitiesSuggestable, filterFullnames, getSuggestionsMatch, clearSuggestions])

  return <Paper className={clsx(className, styles.root)}>
    <InputText
      value={inputValue || null}
      onChange={handleInputChange}
      onSelect={handleAccept}
      onAccept={handleAccept}
      onHighlight={handleHighlight}
      suggestions={keys}
      disableAcceptOnBlur
      autoHighlight={inputValue?.trim() === suggestions[keys[0]]?.input?.trim()}
      ListboxComponent={ListboxSuggestion}
      TextFieldProps={{
        variant: 'outlined',
        placeholder: 'Type your query or keyword here',
        label: error || undefined,
        error: !!error,
        InputLabelProps: { shrink: true },
        size: "medium"
      }}
      InputProps={{
        startAdornment: <SearchIcon className={styles.iconButton} color="action" />,
        endAdornment: <Tooltip title="Search bar syntax help">
          <HelpButton
            IconProps={{fontSize: 'small'}}
            maxWidth="md"
            size="small"
            heading="Search bar help"
            text={`
The search bar provides a fast way to start formulating queries.
Once you start typing a keyword or a query, suggestions for queries
and metainfo names are given based on your search history and the
available data. This search bar supports the following syntaxes:

${formatReadmeList}`}
          />
        </Tooltip>,
        inputRef: inputRef
      }}
      getOptionLabel={option => option}
      filterOptions={(options) => options}
      loading={loading}
      renderOption={(id) => <Suggestion
        suggestion={suggestions[id]}
        onDelete={() => removeSuggestion(suggestions[id].key)}
        tooltip={suggestions[id].type === SuggestionType.Name ? filterData[suggestions[id].input]?.description : undefined}
      />}
    />
  </Paper>
})

SearchBar.propTypes = {
 // List of quantities that the search targets. The suggestions are returned in
 // the same order.
  quantities: PropTypes.arrayOf(PropTypes.shape({
    name: PropTypes.string,
    size: PropTypes.number
  })),
  className: PropTypes.string
}

/**
 * Custom virtualized list component for displaying searchbar suggestions.
 */
const OuterElementContext = createContext({})
export const ListboxSuggestion = React.forwardRef((props, ref) => {
  const { children, ...other } = props
  const itemSize = 48
  const headerSize = 40
  const itemData = React.Children.toArray(children)
  const itemCount = itemData.length

  // Calculate size of child element.
  const getChildSize = (child) => {
    return React.isValidElement(child) && child.type === ListSubheader
      ? headerSize
      : itemSize
  }

  // Calculates the height of the suggestion box
  const getHeight = () => {
    return itemCount > 8
      ? 8 * itemSize
      : itemData.map(getChildSize).reduce((a, b) => a + b, 0)
  }

  const gridRef = useResetCache(itemCount)

  return <div ref={ref}>
    <OuterElementContext.Provider value={other}>
      <List disablePadding>
        <VariableSizeList
          itemData={itemData}
          height={getHeight() + 2 * LISTBOX_PADDING}
          width="100%"
          ref={gridRef}
          outerElementType={OuterElementType}
          innerElementType="ul"
          itemSize={(index) => getChildSize(itemData[index])}
          overscanCount={5}
          itemCount={itemCount}
        >
          {renderRow}
        </VariableSizeList>
      </List>
    </OuterElementContext.Provider>
  </div>
})

ListboxSuggestion.propTypes = {
  children: PropTypes.node
}

const OuterElementType = React.forwardRef((props, ref) => {
  const outerProps = useContext(OuterElementContext)
  return <div ref={ref} {...props} {...outerProps} />
})

export default SearchBar
