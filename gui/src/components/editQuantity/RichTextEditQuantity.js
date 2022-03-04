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
import React, { useCallback, useState } from 'react'
import { Editor } from '@tinymce/tinymce-react'
import PropTypes from 'prop-types'
import { Box, FormControl, FormLabel, makeStyles } from '@material-ui/core'

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
  const initialHeight = 500
  const {quantityDef, section, onChange} = props
  const [text] = useState(section[quantityDef.name])
  const [focus, setFocus] = useState(false)
  const [initialized, setInitialized] = useState(false)

  const handleChange = useCallback((value, editor) => {
    if (onChange) {
      onChange(value === '' ? undefined : value, section, quantityDef)
    }
  }, [onChange, quantityDef, section])

  const handleImageUpload = useCallback((blobInfo, success, failure, progress) => {
    success('data:' + blobInfo.blob().type + ';base64,' + blobInfo.base64())
  }, [])

  const handleEditorInit = useCallback(editor => {
    setInitialized(true)
  }, [setInitialized])

  return (
    <FormControl fullWidth focused={focus} className={focus ? classes.focused : classes.root}>
      <Box marginY={1}>
        <FormLabel>{quantityDef.name}</FormLabel>
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
          initialValue={text}
        />
      </Box>
    </FormControl>
  )
})
RichTextEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

export default RichTextEditQuantity
