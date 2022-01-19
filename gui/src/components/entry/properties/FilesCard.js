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
import {PropertyCard} from './PropertyCard'
import { Box, Grid, Typography } from '@material-ui/core'

const mockupFiles = [
  'https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcRq9I1MSFAtRvqich4XuLYa7mCOdEDSRZTvlaBcFVJKdf8sJvYIgw0OLmRpSQrxAOyI&usqp=CAU',
  'https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcQ1JXPCcqGPtjomnhyugREmlhAUTcaC-zVcw4wQI933ByxZjwjNjheHJVfeqJr5_pbM8rU&usqp=CAU',
  'https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcRMbYxCdklQe9UAD59h3YrSScEn9ChLGLmE7hSdqB4jPFTOX-kAKlg3bVKK5RdH0lFx3dM&usqp=CAU'
]

const FilesCard = React.memo((props) => {
  return <PropertyCard title="Files">
    <Box padding={2}>
      <Grid container spacing={2}>
        {mockupFiles.map((url, index) => (
          <Grid item key={index}>
            <Box display="flex" flexDirection="column" alignItems="center" width={150}>
              <img src={url} width={150} height={150} />
              <Typography variant="caption">
                example-{index}.png
              </Typography>
            </Box>
          </Grid>
        ))}
      </Grid>
    </Box>
  </PropertyCard>
})

FilesCard.propTypes = {}

export default FilesCard
