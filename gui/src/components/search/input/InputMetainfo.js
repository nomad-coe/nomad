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
  useRef,
  useContext,
  createContext
} from 'react'
import { makeStyles } from '@material-ui/core/styles'
import PropTypes from 'prop-types'
import { Tooltip, List, ListItemText, ListSubheader } from '@material-ui/core'
import HelpOutlineIcon from '@material-ui/icons/HelpOutline'
import { getSchemaAbbreviation, getSuggestions } from '../../../utils'
import { useSearchContext } from '../SearchContext'
import { VariableSizeList } from 'react-window'
import { InputText } from './InputText'

/**
 * A metainfo option shown as a suggestion.
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
 * Wrapper around InputText that is specialized in showing metainfo options.
 */
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
  optional,
  TextFieldProps,
  InputProps,
  PaperComponent,
  group,
  loading
}) => {
  // Predefine all option objects, all option paths and also pre-tokenize the
  // options for faster matching.
  const { keys, keysSet, filter } = useMemo(() => {
    const keys = Object.keys(options)
    const keysSet = new Set(keys)
    const { filter } = getSuggestions(keys, 0)
    return { keys, keysSet, filter }
  }, [options])

  // Used to validate the input and raise errors
  const validate = useCallback((value) => {
    const empty = !value || value.length === 0
    if (optional && empty) {
      return {valid: true, error: undefined}
    } else if (empty) {
      return {valid: false, error: 'Please specify a value.'}
    } else if (!(keysSet.has(value))) {
      return {valid: false, error: 'Invalid value for this field.'}
    }
    return {valid: true, error: undefined}
  }, [keysSet, optional])

  // Handles the final acceptance of a value
  const handleAccept = useCallback((key) => {
    const {valid, error} = validate(key)
    if (valid) {
      onAccept && onAccept(key, options[key])
    } else {
      onError && onError(error)
    }
  }, [validate, onError, onAccept, options])

  // Handles the final acceptance of a value
  const handleSelect = useCallback((key) => {
    onSelect && onSelect(key, options[key])
  }, [onSelect, options])

  // Used to filter the shown options based on input
  const filterOptions = useCallback((opt, { inputValue }) => {
    let filtered = filter(inputValue)
    if (group) filtered = filtered.sort((a, b) => options[a].group > options[b].group ? 1 : -1)
    return filtered
  }, [options, filter, group])

  return <InputText
    value={value || null}
    label={label}
    error={error}
    onChange={onChange}
    onSelect={handleSelect}
    onAccept={handleAccept}
    onBlur={onBlur}
    onError={onError}
    suggestions={keys}
    ListboxComponent={ListboxMetainfo}
    TextFieldProps={TextFieldProps}
    PaperComponent={PaperComponent}
    InputProps={InputProps}
    groupBy={group && ((key) => options?.[key]?.group)}
    renderGroup={group && renderGroup}
    getOptionLabel={option => option}
    filterOptions={filterOptions}
    loading={loading}
    renderOption={(id) => <MetainfoOption id={id} options={options} />}
  />
})

InputMetainfo.propTypes = {
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
  loading: PropTypes.bool
}

InputMetainfo.defaultProps = {
  TextFieldProps: {label: "quantity"}
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
    const suggestions = Object.fromEntries(
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
