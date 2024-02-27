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

import React, { useCallback, useContext, useEffect, useMemo, useState } from 'react'
import PropTypes from 'prop-types'
import { isNil } from 'lodash'
import { FilterSubMenu, filterMenuContext } from './FilterMenu'
import {
  Box,
  Button,
  FormControl,
  InputLabel,
  LinearProgress,
  MenuItem,
  Select,
  Typography
} from '@material-ui/core'
import { InputGrid, InputGridItem } from '../input/InputGrid'
import { useApi } from '../../api'
import { useErrors } from '../../errors'
import { useSearchContext } from '../SearchContext'
import { useGlobalMetainfo } from '../../archive/metainfo'
import { NumberEditQuantity } from '../../editQuantity/NumberEditQuantity'
import { StringEditQuantity } from '../../editQuantity/StringEditQuantity'
import { EnumEditQuantity } from '../../editQuantity/EnumEditQuantity'
import { DateTimeEditQuantity, DateEditQuantity } from '../../editQuantity/DateTimeEditQuantity'
import { editQuantityComponents } from '../../editQuantity/EditQuantity'
import { InputMetainfoControlled } from '../../search/input/InputMetainfo'
import { DType, getDatatype, rsplit, parseQuantityName } from '../../../utils'
import { InputTextField } from '../input/InputText'

const types = {
  [DType.Int]: 'Integer',
  [DType.Float]: 'Float',
  [DType.Timestamp]: 'Datetime',
  [DType.String]: 'String',
  [DType.Enum]: 'Enum'
}

const operators = {
  'search': '=',
  'lt': '<',
  'lte': '<=',
  'gte': '=>',
  'gt': '>'
}

const typeOperators = {
  [DType.Int]: Object.keys(operators),
  [DType.Float]: Object.keys(operators),
  [DType.Timestamp]: Object.keys(operators),
  [DType.String]: ['search'],
  [DType.Enum]: ['search'],
  [DType.Boolean]: ['search']
}

const valueKeys = {
  [DType.Int]: 'int_value',
  [DType.Float]: 'float_value',
  [DType.Timestamp]: 'datetime_value',
  [DType.String]: 'str_value',
  [DType.Enum]: 'str_value',
  [DType.Boolean]: 'bool_value'
}

const componentMap = {
  [DType.Int]: NumberEditQuantity,
  [DType.Float]: NumberEditQuantity,
  [DType.Timestamp]: DateEditQuantity,
  [DType.String]: StringEditQuantity,
  [DType.Enum]: EnumEditQuantity,
  [DType.Boolean]: StringEditQuantity
}

const getValueKey = (quantityDef) => {
  const dtype = getDatatype(quantityDef)
  const result = valueKeys[dtype]

  if (!result) {
    throw Error('Unsupported quantity type')
  }
  return result
}

const getOperators = (quantityDef) => {
  const dtype = getDatatype(quantityDef)
  const result = typeOperators[dtype]
  if (!result) {
    throw Error('Unsupported quantity type')
  }
  return result
}

const EditQuantity = React.memo(({quantityDef, value, onChange}) => {
  const [component, componentProps] = useMemo(() => {
    const elnAnnotation = quantityDef.m_annotations?.eln[0] || {}
    let {component, props, ...otherProps} = elnAnnotation
    component = component && editQuantityComponents[component]
    const usesCompatibleComponent = [StringEditQuantity, EnumEditQuantity, NumberEditQuantity, DateTimeEditQuantity].indexOf(component) !== -1

    if (!(component && usesCompatibleComponent)) {
      const type = getDatatype(quantityDef)
      // fall back on StringEditQuantity if type ends up being unknown
      component = type === DType.Unknown ? StringEditQuantity : componentMap[type]
      props = {}
      otherProps = {}
    }

    return [component, {
      ...props,
      ...otherProps,
      label: 'value',
      quantityDef: quantityDef
    }]
  }, [quantityDef])

  return React.createElement(
    component, {
      ...componentProps,
      value: value,
      onChange: onChange
  })
})
EditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.any,
  onChange: PropTypes.func.isRequired
}

const QuantityFilter = React.memo(({quantities, filter, onChange}) => {
  const {id, value} = filter
  const operator = filter.operator || 'search'
  const quantityDef = useMemo(() => {
    return quantities.find(q => q.id === id)?._quantityDef
  }, [id, quantities])

  const options = useMemo(() => {
    return Object.fromEntries(quantities.map(quantity => {
      const {path, schema} = parseQuantityName(quantity.id)
      return [
        quantity.id, {
          key: quantity.id,
          primary: path,
          schema: schema,
          description: quantity?._quantityDef?.description,
          dtype: getDatatype(quantity?._quantityDef)
        }
      ]
    }
      ))
  }, [quantities])

  const handleValueChange = useCallback((value) => {
    onChange({
      operator: 'search',
      ...filter,
      value: value
    })
  }, [filter, onChange])

  const handleIdChange = useCallback((value) => {
    onChange({
      operator: 'search',
      id: value
    })
  }, [onChange])

  const handleOperatorChange = useCallback((e) => {
    onChange({
      ...filter,
      operator: e.target.value
    })
  }, [filter, onChange])

  const availableOperators = quantityDef ? getOperators(quantityDef) : ['search']

  return (<React.Fragment>
    <Box display="flex" flexWrap="wrap" flexDirection="row" alignItems="flex-start" marginTop={1}>
      <Box marginBottom={1} width="100%">
        <InputMetainfoControlled
          options={options}
          value={id}
          onChange={handleIdChange}
          onSelect={handleIdChange}
          // TODO: There should be better error handling here. Errors for now
          // simply clear out the input.
          onError={(error) => error && handleIdChange("")}
          optional
        />
      </Box>
      <Box marginRight={1} width={65}>
        <FormControl fullWidth size="small" variant="filled">
          <InputLabel>op</InputLabel>
          <Select
            value={operator}
            onChange={handleOperatorChange}
            renderValue={value => <b>{operators[value]}</b>}
          >
            {availableOperators.map((operator, index) => (
              <MenuItem key={index} value={operator}>{operators[operator]}</MenuItem>
            ))}
          </Select>
        </FormControl>
      </Box>
      <Box flexGrow={1}>
        {
          quantityDef ? (
            <EditQuantity
              fullWidth
              value={value}
              onChange={handleValueChange}
              quantityDef={quantityDef}
            />
          ) : (
            <InputTextField
              label="value"
              fullWidth
              value={filter.value || ''}
              onChange={e => handleValueChange(e.target.value)}
            />
          )
        }
      </Box>
    </Box>
  </React.Fragment>)
})
QuantityFilter.propTypes = {
  quantities: PropTypes.arrayOf(PropTypes.object).isRequired,
  filter: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

const FilterSubMenuCustomQuantities = React.memo(({
  id,
  ...rest
}) => {
  const metainfo = useGlobalMetainfo()
  const {selected, open} = useContext(filterMenuContext)
  const {useFilterState} = useSearchContext()
  const [loaded, setLoaded] = useState(false)
  const visible = open && id === selected
  const {api} = useApi()
  const {raiseError} = useErrors()
  const [query, setQuery] = useFilterState('custom_quantities')
  const [quantities, setQuantities] = useState(null)

  const [andFilters, setAndFilters] = useState([{}])

  useEffect(() => {
    if (!metainfo || loaded || !visible) {
      return
    }
    const retrieve = async () => {
      const response = await api.post('entries/query', {
        'owner': 'visible',
        'query': {
          'quantities:any': ['data', 'nexus']
        },
        'pagination': {
          'page_size': 0
        },
        'aggregations': {
          'paths': {
            'terms': {
              'quantity': 'search_quantities.id',
              'size': 1000,
              'entries': {
                'size': 1,
                'required': {
                  'include': [
                    'search_quantities.id',
                    'search_quantities.definition'
                  ]
                }
              }
            }
          }
        }
      }, {
        noLoading: true
      })
      const quantities = []
      for (const path of response.aggregations.paths.terms.data) {
        const searchableQuantity = path.entries[0].search_quantities
        const [section_path, quantity_path] = rsplit(searchableQuantity.definition, '.', 1)

        // The definition may not exist in the installation if e.g. the
        // definition is removed, or access rights to YAML files are changed.
        // In this case the value is not shown.
        let quantityDef
        try {
          const sectionDef = await metainfo.resolveDefinition(section_path, false)
          quantityDef = sectionDef?._properties?.[quantity_path]
        } catch (error) {
        }
        if (!quantityDef) continue

        const type = getDatatype(quantityDef)
        let description = types[type] || 'unknown type'
        if (quantityDef.unit) {
          description += ` in ${quantityDef.unit}`
        }
        quantities.push({
          ...searchableQuantity,
          _quantityDef: quantityDef,
          _description: description
        })
      }
      setQuantities(quantities)
      setLoaded(true)
    }
    retrieve().catch(raiseError)
  }, [loaded, visible, api, raiseError, setQuantities, metainfo])

  const handleFilterChange = useCallback((filter, index) => {
    setAndFilters(filters => {
      const newFilters = [...filters]
      newFilters[index] = filter
      return newFilters
    })
  }, [setAndFilters])

  const handleAndClicked = useCallback(() => {
    setAndFilters(filters => [...filters, {}])
  }, [setAndFilters])

  const handleClearClicked = useCallback(() => {
    setAndFilters([{}])
    setQuery([])
  }, [setAndFilters, setQuery])

  const searchEnabled = useMemo(() => {
    return andFilters.every(f => f.id && !isNil(f.value))
  }, [andFilters])

  const handleSearchClicked = useCallback(() => {
    const query = {
      and: andFilters.map(filter => {
        const {id, value, operator} = filter
        const quantityDef = id && quantities.find(q => q.id === id)._quantityDef
        const valueKey = getValueKey(quantityDef)
        const valueKeyWithOperator = operator === 'search' ? valueKey : `${valueKey}:${operator}`
        return {
          search_quantities: {
            id: id,
            [valueKeyWithOperator]: value
          }
        }
      })
    }
    setQuery(query)
  }, [andFilters, setQuery, quantities])

  useEffect(() => {
    const andFilters = query?.and?.map(filter => {
      const searchQuantity = filter.search_quantities
      const valueKey = Object.keys(searchQuantity)[1]
      const valueKeyWithOperator = valueKey.split(':')
      const operator = valueKeyWithOperator.length === 2 ? valueKeyWithOperator[1] : 'search'
      return {
        id: searchQuantity.id,
        value: searchQuantity[Object.keys(searchQuantity)[1]],
        operator: operator
      }
    })
    setAndFilters(andFilters || [{}])
  }, [query, setAndFilters])

  if (!quantities) {
    return (
      <FilterSubMenu id={id} {...rest}>
        <InputGrid>
          <InputGridItem xs={12}>
            <LinearProgress/>
          </InputGridItem>
        </InputGrid>
      </FilterSubMenu>
    )
  }

  if (quantities.length === 0) {
    return (
      <FilterSubMenu id={id} {...rest}>
        <InputGrid>
          <InputGridItem xs={12}>
            <Typography>
              Could not find valid quantities to search for. Make sure that you have access to the data and that the metainfo associated with the quantities can be loaded correctly.
            </Typography>
          </InputGridItem>
        </InputGrid>
      </FilterSubMenu>
    )
  }

  return <FilterSubMenu id={id} {...rest}>
    <InputGrid>
      <InputGridItem xs={12}>
        <Box marginTop={2}>
          {andFilters.map((filter, index) => (
            <Box key={index} marginY={1}>
              <QuantityFilter
                quantities={quantities}
                filter={filter}
                onChange={filter => handleFilterChange(filter, index)}
              />
            </Box>
          ))}
          <Box flexDirection="row" display="flex">
            <Button
              variant="contained" color="primary" onClick={handleAndClicked}
            >And</Button>
            <Box marginLeft={1}>
              <Button
                variant="contained" color="primary" onClick={handleClearClicked}
                disabled={andFilters.length === 1 && Object.keys(andFilters[0]).length === 0}
              >Clear</Button>
            </Box>
          </Box>
        </Box>
      </InputGridItem>
      <InputGridItem xs={12}>
        <Box marginTop={2} display="flex" flexDirection="row" justifyContent="right">
          <Button
            variant="contained" color="primary" onClick={handleSearchClicked}
            disabled={!searchEnabled}
          >Update search</Button>
        </Box>
      </InputGridItem>
    </InputGrid>
  </FilterSubMenu>
})
FilterSubMenuCustomQuantities.propTypes = {
  id: PropTypes.string
}

export default FilterSubMenuCustomQuantities
