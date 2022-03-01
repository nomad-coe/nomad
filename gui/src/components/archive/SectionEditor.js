import React, { useCallback, useState } from 'react'
import PropTypes from 'prop-types'
import { Box, makeStyles, TextField } from '@material-ui/core'
import { useEntryContext } from '../entry/EntryContext'
import _ from 'lodash'
import { AutocompleteEditQuantity, BoolEditQuantity, DateTimeEditQuantity, EnumEditQuantity, NumberEditQuantity, StringEditQuantity } from '../editQuantity/EditQuantity'
import FileEditQuantity from '../editQuantity/FileEditQuantity'
import KeepMaxHeight from '../utils/KeepMaxHeight'
import RichTextEditQuantity from '../editQuantity/RichTextEditQuantity'

const editQuantityComponents = {
  NumberEditQuantity: NumberEditQuantity,
  StringEditQuantity: StringEditQuantity,
  EnumEditQuantity: EnumEditQuantity,
  AutocompleteEditQuantity: AutocompleteEditQuantity,
  BoolEditQuantity: BoolEditQuantity,
  FileEditQuantity: FileEditQuantity,
  DateTimeEditQuantity: DateTimeEditQuantity,
  RichTextEditQuantity: RichTextEditQuantity
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

  return <KeepMaxHeight>
    <div className={classes.root}>
      {showJson ? <JsonEditor data={section} onChange={handleJsonChange} /> : (
        sectionDef._allProperties.map(property => (
          <Box marginBottom={1} key={property.name}>
            <PropertyEditor
              quantityDef={property}
              section={section || {}} onChange={value => handleChange(property, value)}
            />
          </Box>
        ))
      )}
    </div>
  </KeepMaxHeight>
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
