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
import React, {useCallback} from 'react'
import PropTypes from 'prop-types'
import { Box, Grid } from '@material-ui/core'

const ListEditQuantity = React.memo(function ListEditQuantity({value, onChange, component, quantityDef, ...props}) {
  const fixedLength = Number(quantityDef.shape?.[0])
  const hasFixedLength = !isNaN(fixedLength)

  let renderValue
  if (hasFixedLength) {
    renderValue = (value && [...value]) || Array(fixedLength).fill(undefined)
  } else {
    renderValue = (value && [...value]) || []
    renderValue.push(undefined)
  }

  const handleChange = useCallback((item, index) => {
    if (onChange) {
      let newValue = [
        ...(renderValue.slice(0, index)),
        item,
        ...(renderValue.slice(index + 1))
      ]

      if (!hasFixedLength) {
        newValue = newValue.filter(value => value !== undefined && value !== null && value !== '')
      } else {
        for (let i = newValue.length - 1; i >= fixedLength; i--) {
          newValue = newValue.slice(0, i)
        }
      }
      onChange(newValue)
    }
  }, [onChange, renderValue, fixedLength, hasFixedLength])

  return (
    <Box marginY={renderValue.length > 1 ? 2 : undefined}>
      <Grid container direction="column" spacing={1}>
        {renderValue.map((item, index) => (
          <Grid item key={index}>{React.createElement(component, {
            value: item,
            onChange: value => handleChange(value, index),
            quantityDef: quantityDef,
            index: index,
            ...props
          })}</Grid>
        ))}
      </Grid>
    </Box>
  )
})
ListEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.array,
  onChange: PropTypes.func,
  component: PropTypes.elementType.isRequired
}

export default ListEditQuantity
