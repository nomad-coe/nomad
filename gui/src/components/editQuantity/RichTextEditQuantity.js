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
import { useParams } from 'react-router-dom'
import { Editor } from '@tinymce/tinymce-react'
import { useApi } from '../api'
import { useErrors } from '../errors'
import './RichTextEditQuantity.css'

const RichTextEditQuantity = React.memo((props) => {
  const {raiseError} = useErrors()

  const {api} = useApi()
  const { uploadId } = useParams()

  const handleImageUpload = (blobInfo, success, failure, progress) => {
    const formData = new FormData() // eslint-disable-line no-undef
    const fileBlob = blobInfo.blob()
    formData.append('file', fileBlob)

    api.put(`/uploads/${uploadId}/raw/entry_images/`, formData)
      .then(results => console.log(results))
      .catch(raiseError)
      .finally(() => {
        console.log('Uploaded')
      })

    success()
  }
  return (
    <Editor
      init={{
        resize: false,
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
    />
  )
}
)

export default RichTextEditQuantity
