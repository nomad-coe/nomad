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
import React, { useCallback, useState, useRef } from 'react'
import { Editor } from '@tinymce/tinymce-react'
import PropTypes from 'prop-types'
import { Box, FormControl, FormLabel, makeStyles } from '@material-ui/core'
import { getFieldProps } from './StringEditQuantity'
import DOMPurify from 'dompurify'

const useStyle = makeStyles(theme => ({
  root: {
    borderBottom: '1px solid rgba(0, 0, 0, 0.42)',
    marginBottom: 1
  },
  focused: {
    transition: 'border-bottom-color 200ms cubic-bezier(0.4, 0, 0.2, 1) 0ms',
    transitionProperty: 'border-bottom-color',
    transitionDuration: '200ms',
    transitionTimingFunction: 'cubic-bezier(0.4, 0, 0.2, 1)',
    transitionDelay: '0ms',
    borderBottom: `2px solid ${theme.palette.primary.main}`
  }
}))

const RichTextEditQuantity = React.memo((props) => {
  const classes = useStyle()
  const {quantityDef, value, onChange} = props
  const initialHeight = 500
  const {label} = getFieldProps(quantityDef)
  const initialValue = useRef(value)
  const editedValue = useRef(value)
  const [focus, setFocus] = useState(false)
  const [initialized, setInitialized] = useState(false)

  if (editedValue.current !== value) {
    // Means that the value has been changed elsewhere, for example edited in a different tab
    // In this case we need to force an update to the provided value
    editedValue.current = value
    initialValue.current = value
  }

  const handleChange = useCallback((value) => {
    editedValue.current = value
    if (onChange) {
      onChange(value === '' ? undefined : value)
    }
  }, [onChange])

  const handleImageUpload = useCallback((blobInfo, success, failure, progress) => {
    success('data:' + blobInfo.blob().type + ';base64,' + blobInfo.base64())
  }, [])

  const handleEditorInit = useCallback(editor => {
    setInitialized(true)
  }, [setInitialized])

  return (
    <FormControl fullWidth focused={focus} className={focus ? classes.focused : classes.root}>
      <Box marginY={1}>
        <FormLabel>{label}</FormLabel>
      </Box>
      <Box height={initialized ? 'initial' : initialHeight}>
        <Editor
          onInit={(event, editor) => handleEditorInit(editor)}
          init={{
            resize: true,
            height: initialHeight,
            menubar: false,
            plugins: [
              'advlist autolink lists link image charmap print preview anchor',
              'searchreplace visualblocks code',
              'insertdatetime media table paste code help wordcount'
            ],
            toolbar: 'undo redo | formatselect | ' +
            'bold italic backcolor editimage | alignleft aligncenter ' +
            'alignright alignjustify | bullist numlist outdent indent | image table | ' +
            'removeformat',
            skin: 'nomad',
            images_upload_handler: handleImageUpload,
            paste_data_images: true
          }}
          onEditorChange={handleChange}
          onFocus={() => setFocus(true)}
          onBlur={() => setFocus(false)}
          initialValue={DOMPurify.sanitize(initialValue.current || '')}
        />
      </Box>
    </FormControl>
  )
})
RichTextEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  value: PropTypes.string,
  onChange: PropTypes.func
}

export default RichTextEditQuantity
