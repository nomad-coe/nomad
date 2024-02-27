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

import React, {
  useCallback,
  useMemo,
  useEffect,
  useContext,
  createContext,
  useRef
} from 'react'
import { makeStyles } from '@material-ui/core/styles'
import PropTypes from 'prop-types'
import { Tooltip, List, ListItemText, ListSubheader } from '@material-ui/core'
import HelpOutlineIcon from '@material-ui/icons/HelpOutline'
import { getSchemaAbbreviation, getSuggestions, parseJMESPath } from '../../../utils'
import { useSearchContext } from '../SearchContext'
import { VariableSizeList } from 'react-window'
import { InputText } from './InputText'

/**
 * Wrapper around InputText that is specialized in showing metainfo options. The
 * allowed options are controlled.
 */
export const InputMetainfoControlled = React.memo(({
  label,
  value,
  options,
  error,
  validate,
  onChange,
  onSelect,
  onAccept,
  onBlur,
  onError,
  optional,
  TextFieldProps,
  InputProps,
  PaperComponent,
  group,
  loading,
  disableValidation,
  disableValidateOnSelect
}) => {
  // Predefine all option objects, all option paths and also pre-tokenize the
  // options for faster matching.
  const { keys, keysSet, filter } = useMemo(() => {
    const keys = Object.keys(options)
    const keysSet = new Set(keys)
    const { filter } = getSuggestions(keys, 0)
    return { keys, keysSet, filter }
  }, [options])

  const textProps = useMemo(() => {
    return {
      label,
      ...TextFieldProps
    }
  }, [label, TextFieldProps])

  // Used to validate the input and raise errors
  const validateFinal = useCallback((value) => {
    if (disableValidation) {
      return {valid: true, error: undefined}
    }
    const empty = !value || value.length === 0
    if (optional && empty) {
      return {valid: true, error: undefined}
    } else if (empty) {
      return {valid: false, error: 'Please specify a value'}
    }
    if (validate) {
      return validate(value)
    } else if (!(keysSet.has(value))) {
      return {valid: false, error: 'Invalid value for this field'}
    }
    return {valid: true, error: undefined}
  }, [validate, keysSet, optional, disableValidation])

  // Handles the selectance of a suggested value
  const handleSelect = useCallback((key) => {
    onSelect?.(key, options[key])
  }, [onSelect, options])

  // Used to filter the shown options based on input
  const filterOptions = useCallback((opt, { inputValue }) => {
    let filtered = filter(inputValue)
    if (group) filtered = filtered.sort((a, b) => options[a].group > options[b].group ? 1 : -1)
    return filtered
  }, [options, filter, group])

  return <InputText
    value={value || null}
    error={error}
    onChange={onChange}
    onSelect={handleSelect}
    onAccept={onAccept}
    onBlur={onBlur}
    onError={onError}
    suggestions={keys}
    ListboxComponent={ListboxMetainfo}
    TextFieldProps={textProps}
    PaperComponent={PaperComponent}
    InputProps={InputProps}
    groupBy={group && ((key) => options?.[key]?.group)}
    renderGroup={group && renderGroup}
    getOptionLabel={option => option}
    filterOptions={filterOptions}
    loading={loading}
    renderOption={(id) => <MetainfoOption id={id} options={options} />}
    validate={validateFinal}
    disableValidateOnSelect={disableValidateOnSelect}
  />
})

InputMetainfoControlled.propTypes = {
  label: PropTypes.string,
  value: PropTypes.string,
  options: PropTypes.objectOf(PropTypes.shape({
    key: PropTypes.string.isRequired, // The value used for matching
    definition: PropTypes.object, // The definition for the metainfo.
    primary: PropTypes.string, // The value shown as primary text, defaults to key.
    secondary: PropTypes.string, // Optional secondary text. Defaults to data type + schema read from definition.
    description: PropTypes.string, // Optional description shown as hover. Defaults to desription read from definition.
    schema: PropTypes.string, // Schema name. Defaults to schema read from searchQuantities.
    dtype: PropTypes.string, // Data type. Defaults to data type read from definition.
    group: PropTypes.string // Optional group information
  })),
  error: PropTypes.string,
  validate: PropTypes.func, // Optional custom validation function
  onChange: PropTypes.func,
  onSelect: PropTypes.func,
  onAccept: PropTypes.func,
  onError: PropTypes.func,
  onBlur: PropTypes.func,
  optional: PropTypes.bool, // Set to true if field can be empty
  group: PropTypes.bool, // Set to true if group should be shown
  TextFieldProps: PropTypes.object,
  InputProps: PropTypes.object,
  PaperComponent: PropTypes.any,
  loading: PropTypes.bool,
  disableValidation: PropTypes.bool,
  disableValidateOnSelect: PropTypes.bool
}

InputMetainfoControlled.defaultProps = {
  label: "quantity"
}

/**
 * Wrapper around InputText that is specialized in showing metainfo options. The
 * allowed options are uncontrolled: they are specified by the current
 * SearchContext.
 */
export const InputMetainfo = React.memo(({
  dtypes,
  dtypesRepeatable,
  disableNonAggregatable,
  ...rest
}) => {
  const { filterData } = useSearchContext()

  // Fetch the available metainfo names and create options that are compatible
  // with InputMetainfo.
  const options = useMemo(() => {
    return getMetainfoOptions(filterData, dtypes, dtypesRepeatable, disableNonAggregatable)
  }, [filterData, dtypes, dtypesRepeatable, disableNonAggregatable])

  return <InputMetainfoControlled
    options={options}
    {...rest}
  />
})

InputMetainfo.propTypes = {
  /* List of allowed data types for non-repeatable quantities. */
  dtypes: PropTypes.object,
  /* List of allowed data types for repeatable quantities. */
  dtypesRepeatable: PropTypes.object,
  /* Whether non-aggregatable values are excluded */
  disableNonAggregatable: PropTypes.bool
}

/**
 * Wrapper around InputText which accepts JMESPath syntax for quantities. The
 * allowed quantities are uncontrolled: they are specified by the current
 * SearchContext.
 */
export const InputJMESPath = React.memo(React.forwardRef(({
  dtypes,
  dtypesRepeatable,
  disableNonAggregatable,
  ...rest
}, ref) => {
  const { filterData } = useSearchContext()

  // Fetch the available metainfo names and create options that are compatible
  // with InputMetainfo.
  const [options, keysSet] = useMemo(() => {
    const options = getMetainfoOptions(filterData, dtypes, dtypesRepeatable, disableNonAggregatable)
    const keysSet = new Set(Object.keys(options))
    return [options, keysSet]
  }, [filterData, dtypes, dtypesRepeatable, disableNonAggregatable])

  // Used to validate the JMESPath input
  const validate = useCallback((value) => {
    const {quantity, path, extras, error: errorParse, schema} = parseJMESPath(value)
    if (errorParse) {
      return {valid: false, error: 'Invalid JMESPath query, please check your syntax.'}
    }
    if (!(keysSet.has(quantity))) {
      return {valid: false, error: `The quantity "${quantity}" is not available.`}
    }
    for (const extra of extras) {
      if (!filterData[extra]) {
        return {valid: false, error: `The quantity "${extra}" is not available.`}
      }
    }
    if (filterData[quantity].repeats_section && quantity === path + schema) {
      return {valid: false, error: `The quantity "${quantity}" is contained in at least one repeatable section. Please use JMESPath syntax to select one or more target sections.`}
    }
    return {valid: true, error: undefined}
  }, [keysSet, filterData])

  return <InputMetainfoControlled
    options={options}
    validate={validate}
    disableValidateOnSelect
    ref={ref}
    {...rest}
  />
}))

InputJMESPath.propTypes = {
  /* List of allowed data types for non-repeatable quantities. */
  dtypes: PropTypes.object,
  /* List of allowed data types for repeatable quantities. */
  dtypesRepeatable: PropTypes.object,
  /* Whether non-aggregatable values are excluded */
  disableNonAggregatable: PropTypes.bool
}

/**
 * Custom virtualized list component for displaying metainfo values.
 */
export const ListboxMetainfo = React.forwardRef((props, ref) => {
  const { children, ...other } = props
  const itemSize = 64
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

ListboxMetainfo.propTypes = {
  children: PropTypes.node
}

export const LISTBOX_PADDING = 8
const OuterElementContext = createContext({})

export const OuterElementType = React.forwardRef((props, ref) => {
  const outerProps = useContext(OuterElementContext)
  return <div ref={ref} {...props} {...outerProps} />
})

export const renderGroup = (params) => [
  <ListSubheader key={params.key} component="div">
    {`${params.group} suggestions`}
  </ListSubheader>,
  params.children
]

export function useResetCache(data) {
  const ref = useRef(null)
  useEffect(() => {
    if (ref.current != null) {
      ref.current.resetAfterIndex(0, true)
    }
  }, [data])
  return ref
}

export function renderRow({ data, index, style }) {
  return React.cloneElement(data[index], {
    style: {
      ...style,
      top: style.top + LISTBOX_PADDING
    }
  })
}

/* A metainfo option shown as a suggestion.
 */
export const useMetainfoOptionStyles = makeStyles(theme => ({
  optionText: {
    flexGrow: 1,
    overflowX: 'scroll',
    '&::-webkit-scrollbar': {
      display: 'none'
    },
    '-ms-overflow-style': 'none',
    scrollbarWidth: 'none'
  },
  noWrap: {
    whiteSpace: 'nowrap'
  },
  option: {
    width: '100%',
    display: 'flex',
    alignItems: 'stretch',
    // The description icon is hidden until the item is hovered. It is not
    // removed from the document with "display: none" in order for the hover to
    // not change the layout which may cause other elements to shift around.
    '& .description': {
      visibility: "hidden",
      marginRight: theme.spacing(-0.5),
      display: 'flex',
      width: theme.spacing(5),
      marginLeft: theme.spacing(1),
      alignItems: 'center',
      justifyContent: 'center'
    },
    '&:hover .description': {
      visibility: "visible"
    }
  }
}))
export const MetainfoOption = ({id, options}) => {
  const styles = useMetainfoOptionStyles()
  const option = options[id]
  const primary = option.primary || option.definition?.quantity || option.key
  const dtype = option.dtype || option.definition?.dtype
  const schema = getSchemaAbbreviation(option.schema || option.definition?.schema)
  const secondary = option.secondary || ((dtype || schema)
    ? `${dtype} ${schema ? `| ${schema}` : ''}`
    : null)
  const description = option.description || option.definition?.description

  return <div className={styles.option}>
    <ListItemText
      primary={primary}
      secondary={secondary}
      className={styles.optionText}
      primaryTypographyProps={{className: styles.noWrap}}
      secondaryTypographyProps={{className: styles.noWrap}}
    />
    {description &&
      <Tooltip title={description || ''}>
        <div className="description">
          <HelpOutlineIcon fontSize="small" color="action"/>
        </div>
      </Tooltip>
    }
  </div>
}

MetainfoOption.propTypes = {
  id: PropTypes.string,
  options: PropTypes.object
}

/**
 * Used to filter a set of filterData according to the given filtering options.
 *
 * @param {*} filterData The filterData to filter
 * @param {*} dtypes Included data types as a set
 * @param {*} dtypesRepeatable  Included repeatabel data types as a set
 * @param {*} disableNonAggregatable  Whether to filter out non-aggregatable filters
 * @returns Object containing the filtered filters
 */
function getMetainfoOptions(filterData, dtypes, dtypesRepeatable, disableNonAggregatable) {
  return Object.fromEntries(
    Object.entries(filterData)
      .filter(([key, data]) => {
        if (disableNonAggregatable && !data.aggregatable) return false
        const dtype = data?.dtype
        return data?.repeats
          ? dtypesRepeatable?.has(dtype)
          : dtypes?.has(dtype)
      })
      .map(([key, data]) => [key, {
        key: key,
        definition: data
    }])
  )
}
