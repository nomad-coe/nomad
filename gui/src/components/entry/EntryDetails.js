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
import React from 'react'
import PropTypes from 'prop-types'
import Quantity from '../Quantity'
import { Formula } from './properties/MaterialCard'

export const MethodDetails = React.memo(({data}) => {
  const methodQuantities = []
  const addMethodQuantities = (obj, parentKey) => {
    const children = {}
    Object.keys(obj).forEach(key => {
      const value = obj[key]
      if (Array.isArray(value) || typeof value === 'string') {
        if (value.length > 0) {
          methodQuantities.push({
            quantity: `${parentKey}.${key}`,
            label: key.replace(/_/g, ' ')
          })
        }
      } else if (value instanceof Object) {
        children[key] = value
      }
    })
    Object.keys(children).forEach(key => addMethodQuantities(children[key], `${parentKey}.${key}`))
  }
  if (data?.results?.method) {
    addMethodQuantities(data.results.method, 'results.method')
  }

  return <Quantity flex>
    {methodQuantities.map(({...quantityProps}) => (
      <Quantity
        key={quantityProps.quantity}
        {...quantityProps}
        noWrap
        data={data}
        hideIfUnavailable
      />
    ))}
  </Quantity>
})

MethodDetails.propTypes = {
  data: PropTypes.object
}

const EntryDetails = React.memo(({data}) => {
  return <>
    <Quantity flex>
      <Formula data={data} />
      <Quantity quantity="results.material.material_name" data={data} label="material name" />
    </Quantity>
    <MethodDetails data={data} />
  </>
})

EntryDetails.propTypes = {
  data: PropTypes.object
}

export default EntryDetails
