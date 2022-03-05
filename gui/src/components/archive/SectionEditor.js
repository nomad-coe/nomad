import React, { useCallback, useMemo, useRef, useState } from 'react'
import PropTypes from 'prop-types'
import { Box, makeStyles, TextField } from '@material-ui/core'
import { useEntryContext } from '../entry/EntryContext'
import _ from 'lodash'
import {DateTimeEditQuantity} from '../editQuantity/DateTimeEditQuantity'
import {StringEditQuantity} from '../editQuantity/StringEditQuantity'
import {NumberEditQuantity} from '../editQuantity/NumberEditQuantity'
import {EnumEditQuantity} from '../editQuantity/EnumEditQuantity'
import {AutocompleteEditQuantity} from '../editQuantity/AutocompleteEditQuantity'
import {BoolEditQuantity} from '../editQuantity/BoolEditQuantity'
import FileEditQuantity from '../editQuantity/FileEditQuantity'
import RichTextEditQuantity from '../editQuantity/RichTextEditQuantity'
import ReferenceEditQuantity from '../editQuantity/ReferenceEditQuantity'
import { QuantityMDef } from './metainfo'

const editQuantityComponents = {
  NumberEditQuantity: NumberEditQuantity,
  StringEditQuantity: StringEditQuantity,
  EnumEditQuantity: EnumEditQuantity,
  SelectEnumEditQuantity: EnumEditQuantity,
  RadioEnumEditQuantity: EnumEditQuantity,
  AutocompleteEditQuantity: AutocompleteEditQuantity,
  BoolEditQuantity: BoolEditQuantity,
  FileEditQuantity: FileEditQuantity,
  DateTimeEditQuantity: DateTimeEditQuantity,
  RichTextEditQuantity: RichTextEditQuantity,
  ReferenceEditQuantity: ReferenceEditQuantity
}

const JsonEditor = React.memo(function JsonEditor({data, onChange}) {
  const [json, setJson] = useState(JSON.stringify(data, null, 2))
  const [error, setError] = useState(null)

  const handleChange = useCallback((event) => {
    const value = event.target.value
    setJson(value)
    try {
      const data = JSON.parse(value)
      if (onChange) {
        onChange(data)
      }
      setError(null)
    } catch (e) {
      setError('This is not JSON: ' + e)
    }
  }, [onChange, setJson])

  return (
    <TextField
      fullWidth label="JSON" error={!!error}
      helperText={error}
      variant="filled" multiline maxRows={20}
      value={json} onChange={handleChange}
    />
  )
})
JsonEditor.propTypes = {
  data: PropTypes.object.isRequired,
  onChange: PropTypes.func
}

const PropertyEditor = React.memo(function PropertyEditor({quantityDef, section, onChange}) {
  const handleChange = useCallback((value) => {
    if (onChange) {
      onChange(value)
    }
  }, [onChange])
  const editAnnotations = quantityDef.m_annotations?.eln || []
  const editAnnotation = editAnnotations[0]
  const componentName = editAnnotation?.component
  const component = componentName && editQuantityComponents[componentName]
  if (!component) {
    return ''
  }
  const props = {
    section: section,
    quantityDef: quantityDef,
    onChange: handleChange,
    ...(editAnnotation?.props || {})
  }
  return React.createElement(component, props)
})
PropertyEditor.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func
}

const useSectionEditorStyles = makeStyles(theme => ({
  root: {
    minWidth: 600
  }
}))
const SectionEditor = React.memo(function SectionEditor({sectionDef, section, onChange, showJson}) {
  const classes = useSectionEditorStyles()
  const {handleArchiveChanged} = useEntryContext()
  const rootRef = useRef()

  const handleChange = useCallback((property, value) => {
    section[property.name] = value
    if (onChange) {
      onChange(section)
    }
    handleArchiveChanged()
  }, [section, onChange, handleArchiveChanged])

  const handleJsonChange = useCallback((data) => {
    _.extend(section, data)
    if (onChange) {
      onChange(section)
    }
    handleArchiveChanged()
  }, [handleArchiveChanged, onChange, section])

  const jsonData = useMemo(() => {
    if (!showJson) {
      return null
    }
    const jsonData = {}
    sectionDef._allProperties
      .filter(property => property.m_def === QuantityMDef && property.m_annotations?.eln)
      .forEach(property => {
        jsonData[property.name] = section[property.name]
      })
    return jsonData
  }, [showJson, section, sectionDef])

  return (
    <div className={classes.root} ref={rootRef}>
      {showJson
        ? (
          <Box height={rootRef.current?.clientHeight} marginY={1}>
            <JsonEditor data={jsonData} onChange={handleJsonChange} />
          </Box>
        ) : (
          sectionDef._allProperties.map(property => (
            <Box marginBottom={1} key={property.name}>
              <PropertyEditor
                quantityDef={property}
                section={section || {}} onChange={value => handleChange(property, value)}
              />
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
