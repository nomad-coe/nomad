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
import { CircularProgress, Grid, makeStyles, Paper, Slide, Typography } from '@material-ui/core'
import Page from '../Page'
import { useUploadPageContext } from './UploadPageContext'

const useStyles = makeStyles(theme => ({
  root: {
    position: 'absolute',
    left: 0,
    right: 0,
    marginTop: -20,
    margin: -10,
    padding: 10,
    paddingTop: 20,
    zIndex: 1000
  }
}))

const UploadProcessingStatus = React.memo(function ProcessingStatus() {
  const {isProcessing, upload} = useUploadPageContext()
  const classes = useStyles()

  if (!upload) {
    return null
  }

  return (
    <Slide direction="down" in={isProcessing} mountOnEnter unmountOnExit>
      <Paper className={classes.root}>
        <Page limitedWidth>
          <Grid container spacing={2} alignItems="center">
            <Grid item>
              <CircularProgress />
            </Grid>
            <Grid item style={{flexGrow: 1}}>
              <Typography>Upload is processing ...</Typography>
              <Typography>{upload.last_status_message}</Typography>
            </Grid>
          </Grid>
        </Page>
      </Paper>
    </Slide>
  )
})

export default UploadProcessingStatus
