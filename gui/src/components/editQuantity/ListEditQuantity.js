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
import React, {useCallback, useState} from 'react'
import PropTypes from 'prop-types'
import {Box, makeStyles} from '@material-ui/core'
import {configState, FoldableList} from '../archive/ArchiveBrowser'
import grey from "@material-ui/core/colors/grey"
import {getDisplayLabel} from "../../utils"
import {useRecoilValue} from "recoil"

const useListEditQuantityStyles = makeStyles(theme => ({
  root: {
    margin: `0 -${theme.spacing(1)}px`,
    padding: `0 ${theme.spacing(1)}px`,
    width: '100%',
    display: 'flex',
    flexDirection: 'row',
    wrap: 'nowrap',
    alignItems: 'center'
  },
  title: {
    flexGrow: 1,
    color: grey[600],
    textDecoration: 'none',
    margin: `0 -${theme.spacing(1)}px`,
    whiteSpace: 'nowrap',
    display: 'flex',
    fontWeight: 'bold'
  },
  selected: {
    whiteSpace: 'nowrap'
  },
  unSelected: {
    '&:hover': {
      backgroundColor: grey[300]
    }
  },
  actions: {}
}))

const ListEditQuantity = React.memo(function ListEditQuantity({value, onChange, component, quantityDef, ...props}) {
  const classes = useListEditQuantityStyles()
  const [open, setOpen] = useState(true)
  const fixedLength = Number(quantityDef.shape?.[0])
  const hasFixedLength = !isNaN(fixedLength)
  const config = useRecoilValue(configState)
  const label = getDisplayLabel(quantityDef, true, config?.showMeta)

  let renderValue
  if (hasFixedLength) {
    renderValue = (value && [...value]) || Array(fixedLength).fill(undefined)
  } else {
    renderValue = (value && [...value]) || []
    renderValue.push(undefined)
  }

  const handleChange = useCallback((item, index) => {
    if (onChange) {
      let newValue = Array.isArray(item) ? [
        ...(renderValue.slice(0, index)),
        ...item,
        ...(renderValue.slice(index + 1))
      ] : [
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

  const handleClick = useCallback(() => {
    setOpen(open => !open)
  }, [])

  return <Box marginTop={1} marginBottom={1}>
    <FoldableList
      label={label}
      className={classes}
      open={open}
      onClick={handleClick}
      nTop={10}
      nBottom={0}
      pageSize={10}
    >
      {renderValue.map((item, index) => (
        <Box marginTop={1} marginLeft={2} marginRight={2} key={index}>{React.createElement(component, {
          value: item,
          onChange: value => handleChange(value, index),
          quantityDef: quantityDef,
          index: index,
          ...props
        })}</Box>
      ))}
    </FoldableList>
  </Box>
})
ListEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.array,
  onChange: PropTypes.func,
  component: PropTypes.elementType.isRequired
}

export default ListEditQuantity
