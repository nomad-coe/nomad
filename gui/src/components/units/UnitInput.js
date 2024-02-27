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

import React, {useCallback, useEffect, useMemo, useRef, useContext, createContext} from 'react'
import PropTypes from 'prop-types'
import {isNil} from 'lodash'
import {getSuggestions} from '../../utils'
import {unitMap} from './UnitContext'
import {parseQuantity} from './Quantity'
import {List, ListItemText, ListSubheader, makeStyles} from '@material-ui/core'
import {VariableSizeList} from 'react-window'
import {InputText} from '../search/input/InputText'

/**
 * Wrapper around InputText that is specialized in showing unit options.
 */
export const useInputStyles = makeStyles(theme => ({
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
    alignItems: 'stretch'
  }
}))
export const UnitInput = React.memo(({value, error, onChange, onAccept, onSelect, onError, dimension, options, disabled, label, disableGroup, optional}) => {
  const styles = useInputStyles()

  // Predefine all option objects, all option paths and also pre-tokenize the
  // options for faster matching.
  const {keys, filter, finalOptions} = useMemo(() => {
    const finalOptions = {}
    Object.entries(unitMap)
      .filter(([key, unit]) => dimension ? unit.dimension === dimension : true)
      .forEach(([key, unit]) => {
        const unitAbbreviation = unit.abbreviation ? ` (${unit.abbreviation})` : ''
        finalOptions[key] = {
          key: key,
          primary: `${unit.label}${unitAbbreviation}`,
          secondary: unit.aliases?.splice(1).join(', '),
          dimension: unit.dimension,
          unit: unit
        }
      })
    const keys = Object.keys(finalOptions)
    const {filter} = getSuggestions(keys, 0)
    return {keys, filter, finalOptions}
  }, [dimension])

  const validate = useCallback((value) => {
    if (isNil(value) || value?.trim?.() === '') {
      if (optional) {
        return {valid: true}
      } else {
        return {valid: false, error: 'Please specify a value'}
      }
    }
    const {error, unit} = parseQuantity(value, dimension, false, true)
    return {valid: !error, error, data: unit}
  }, [optional, dimension])

  // Revalidate input when dimension changes
  useEffect(() => {
    if (!value) return
    const {error} = validate(value)
    if (error) {
      onError?.(error)
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [dimension])

  // Used to filter the shown options based on input
  const filterOptions = useCallback((opt, { inputValue }) => {
    let filtered = filter(inputValue)
    if (!disableGroup) filtered = filtered.sort((a, b) => options[a].group > options[b].group ? 1 : -1)
    return filtered
  }, [disableGroup, filter, options])

  return <InputText
    value={value}
    TextFieldProps={{label, disabled}}
    error={error}
    onChange={onChange}
    onSelect={onSelect}
    onAccept={onAccept}
    onError={onError}
    validate={validate}
    disableClearable
    suggestAllOnFocus
    showOpenSuggestions
    suggestions={keys}
    ListboxComponent={ListboxUnit}
    groupBy={disableGroup ? undefined : (key) => finalOptions?.[key]?.dimension}
    renderGroup={disableGroup ? undefined : renderGroup}
    getOptionLabel={option => option}
    filterOptions={filterOptions}
    renderOption={(key) => {
      const option = finalOptions[key]
      return <div className={styles.option}>
        <ListItemText
          primary={option.primary || option.key}
          secondary={option.secondary}
          className={styles.optionText}
          primaryTypographyProps={{className: styles.noWrap}}
          secondaryTypographyProps={{className: styles.noWrap}}
        />
      </div>
    }}
  />
})
UnitInput.propTypes = {
  value: PropTypes.string,
  error: PropTypes.string,
  label: PropTypes.string,
  dimension: PropTypes.string,
  options: PropTypes.object,
  onChange: PropTypes.func,
  onAccept: PropTypes.func,
  onSelect: PropTypes.func,
  onBlur: PropTypes.func,
  onError: PropTypes.func,
  disabled: PropTypes.bool,
  disableGroup: PropTypes.bool,
  optional: PropTypes.bool
}

export default UnitInput

/**
 * Custom virtualized list component for displaying unit values.
 */
const ListboxUnit = React.forwardRef((props, ref) => {
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

ListboxUnit.propTypes = {
  children: PropTypes.node
}

const LISTBOX_PADDING = 8
const OuterElementContext = createContext({})

const OuterElementType = React.forwardRef((props, ref) => {
  const outerProps = useContext(OuterElementContext)
  return <div ref={ref} {...props} {...outerProps} />
})

const renderGroup = (params) => [
  <ListSubheader key={params.key} component="div">
    {params.group}
  </ListSubheader>,
  params.children
]

function useResetCache(data) {
  const ref = useRef(null)
  useEffect(() => {
    if (ref.current != null) {
      ref.current.resetAfterIndex(0, true)
    }
  }, [data])
  return ref
}

function renderRow({ data, index, style }) {
  return React.cloneElement(data[index], {
    style: {
      ...style,
      top: style.top + LISTBOX_PADDING
    }
  })
}
