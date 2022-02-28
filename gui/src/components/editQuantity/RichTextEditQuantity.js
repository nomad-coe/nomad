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

const RichTextEditQuantity = React.memo((props) => {
  const {quantityDef, section, onChange} = props
  const [text] = useState(section[quantityDef.name])

  const handleChange = useCallback((value, editor) => {
    if (onChange) {
      onChange(value, section, quantityDef)
    }
  }, [onChange, quantityDef, section])

  const handleImageUpload = (blobInfo, success, failure, progress) => {
    success('data:' + blobInfo.blob().type + ';base64,' + blobInfo.base64())
  }
  return (
    <Editor
      init={{
        resize: true,
        height: 500,
        menubar: false,
        plugins: [
          'advlist autolink lists link image charmap print preview anchor',
          'searchreplace visualblocks code',
          'insertdatetime media table paste code help wordcount'
        ],
        toolbar: 'undo redo | formatselect | ' +
        'bold italic backcolor editimage | alignleft aligncenter ' +
        'alignright alignjustify | bullist numlist outdent indent | ' +
        'removeformat',
        skin: 'nomad',
        images_upload_handler: handleImageUpload,
        paste_data_images: true
      }}
      onEditorChange={handleChange}
      initialValue={text}
    />
  )
})
RichTextEditQuantity.propTypes = {
  quantityDef: PropTypes.object.isRequired,
  section: PropTypes.object.isRequired,
  onChange: PropTypes.func.isRequired
}

export default RichTextEditQuantity
