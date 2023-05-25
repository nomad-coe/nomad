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
import { Typography, Box, DialogTitle, DialogContent, Dialog, LinearProgress} from '@material-ui/core'

function LinearProgressWithLabel(props) {
  return (
    <Box display="flex" alignItems="center">
      <Box width="100%" mr={1}>
        <LinearProgress variant="determinate" {...props} />
      </Box>
      <Typography variant="body2" color="textSecondary">{`${Math.round(
        props.value
      )}%`}</Typography>
    </Box>
  )
}
LinearProgressWithLabel.propTypes = {
  value: PropTypes.number.isRequired
}

export default function UploadProgressDialog({uploading}) {
    if (uploading || uploading === 0) {
      return (<Dialog open>
        <DialogTitle>Uploading file ...</DialogTitle>
        <DialogContent>
          <Box width={300}>
            <LinearProgressWithLabel value={uploading} />
          </Box>
        </DialogContent>
      </Dialog>)
  }
  return null
}

UploadProgressDialog.propTypes = {
  uploading: PropTypes.number.isRequired
}
