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

import React, { useCallback, useMemo, useRef } from 'react'
import PropTypes from 'prop-types'
import { Box, makeStyles } from '@material-ui/core'
import { useEntryStore } from '../entry/EntryContext'
import { extend } from 'lodash'
import ListEditQuantity from '../editQuantity/ListEditQuantity'
import InputConfig from '../search/input/InputConfig'
import { editQuantityComponents } from '../editQuantity/EditQuantity'
import { QuantityMDef } from './metainfo'
import {getAllVisibleProperties} from './ArchiveBrowser'
import {QuantityRow, QuantityTable} from '../Quantity'
import {PropertyPreview} from '../entry/properties/SectionCard'

const PropertyEditor = React.memo(function PropertyEditor({quantityDef, value, onChange}) {
  const editAnnotations = quantityDef.m_annotations?.eln || []
  const editAnnotation = editAnnotations[0] || {}
  const {component, props, ...moreProps} = editAnnotation
  const editComponent = component && editQuantityComponents[component]
  if (!editComponent) {
    return null
  }
  const editComponentProps = {
    quantityDef: quantityDef,
    value: value === undefined ? quantityDef.default : value,
    onChange: onChange,
    ...moreProps,
    ...(props || {})
  }

  const shape = quantityDef.shape || []
  if (shape.length === 0) {
    return <Box data-testid={"editable-quantity-editor"}>
      {React.createElement(editComponent, editComponentProps)}
    </Box>
  } else if (shape.length === 1) {
    return <Box data-testid={"editable-quantity-editor"}>
      <ListEditQuantity
        component={editComponent}
        {...editComponentProps}
      />
    </Box>
  } else {
    console.log('Unsupported quantity shape ', shape)
    return null
  }
})
PropertyEditor.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.any,
  onChange: PropTypes.func.isRequired
}

const useSectionEditorStyles = makeStyles(theme => ({
  root: {
    minWidth: 600
  },
  quantityTable: {
    marginBottom: 5
  }
}))
const SectionEditor = React.memo(function SectionEditor({sectionDef, section, onChange, showJson}) {
  const classes = useSectionEditorStyles()
  const {handleArchiveChanged} = useEntryStore() || {}
  const rootRef = useRef()

  const handleChange = useCallback((property, value) => {
    if (section[property.name] === value) {
      return
    }
    section[property.name] = value
    if (onChange) {
      onChange(section)
    }
    handleArchiveChanged()
  }, [section, onChange, handleArchiveChanged])

  const handleJsonChange = useCallback((data) => {
    extend(section, data)
    if (onChange) {
      onChange(section)
    }
    handleArchiveChanged()
  }, [handleArchiveChanged, onChange, section])

  const allVisibleProperties = useMemo(() => getAllVisibleProperties(sectionDef), [sectionDef])
  const allVisibleQuantities = useMemo(() => allVisibleProperties.filter(property => property.m_def === QuantityMDef && property.m_annotations?.eln), [allVisibleProperties])

  const jsonData = useMemo(() => {
    if (!showJson) {
      return null
    }
    const jsonData = {}
    allVisibleQuantities
      .filter(property => {
        // TODO this is just a hack to avoid large values, e.g. rich text with images
        const value = section[property.name]
        return !value || typeof value !== 'string' || value.length <= 1e3
      })
      .forEach(property => {
        jsonData[property.name] = section[property.name]
      })
    return jsonData
  }, [showJson, allVisibleQuantities, section])

  return (
    <div className={classes.root} ref={rootRef}>
      {showJson
        ? (
          <Box height={rootRef.current?.clientHeight} marginY={1}>
            <InputConfig data={jsonData} onChange={handleJsonChange} />
          </Box>
        ) : (
          allVisibleQuantities.map(quantity => (
            quantity._isEditable
              ? <Box marginBottom={1} key={quantity.name} data-testid={"visible-quantity"}>
                <PropertyEditor
                  quantityDef={quantity}
                  value={section?.[quantity.name]} onChange={value => handleChange(quantity, value)}
                />
              </Box>
              : <Box key={quantity.name} data-testid={"visible-quantity"}>
                <QuantityTable className={classes.quantityTable} data={section}>
                  <QuantityRow>
                    <PropertyPreview quantityDef={quantity} section={section}/>
                  </QuantityRow>
                </QuantityTable>
              </Box>
          ))
        )
      }
    </div>
  )
})
SectionEditor.propTypes = {
  sectionDef: PropTypes.object.isRequired,
  section: PropTypes.object,
  onChange: PropTypes.func,
  showJson: PropTypes.bool,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
}

export default SectionEditor
