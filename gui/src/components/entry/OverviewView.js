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
import React, { useState, useEffect, useMemo } from 'react'
import PropTypes from 'prop-types'
import { useApi } from '../api'
import { useErrors } from '../errors'
import { Typography, makeStyles, Box, Grid, Divider } from '@material-ui/core'
import Quantity from '../Quantity'
import ElectronicPropertiesCard from '../entry/properties/ElectronicPropertiesCard'
import MaterialCard from '../entry/properties/MaterialCard'
import VibrationalPropertiesCard from '../entry/properties/VibrationalPropertiesCard'
import MechanicalPropertiesCard from '../entry/properties/MechanicalPropertiesCard'
import GeometryOptimizationCard from '../entry/properties/GeometryOptimizationCard'
import SpectroscopyCard from './properties/SpectroscopyCard'
import { MethodMetadata } from './EntryDetails'
import Page from '../Page'
import { SourceApiCall, SourceApiDialogButton, SourceDialogDivider } from '../buttons/SourceDialogButton'

function MetadataSection({title, children}) {
  return <Box marginTop={2} marginBottom={2}>
    {title && <Typography component="div">
      <Box fontSize="h6.fontSize" marginBottom={1}>
        {title}
      </Box>
    </Typography>}
    {children}
  </Box>
}

MetadataSection.propTypes = {
  title: PropTypes.string,
  children: PropTypes.any
}

const useStyles = makeStyles(theme => ({
  root: {
    marginBottom: theme.spacing(4)
  },
  leftColumn: {
    maxWidth: '32%',
    flexBasis: '32%',
    flexGrow: 0,
    paddingRight: theme.spacing(3)
  },
  rightColumn: {
    maxWidth: '67.99%',
    flexBasis: '67.99%',
    flexGrow: 0,
    '& > div': {
      marginBottom: theme.spacing(2)
    }
  },
  divider: {
    marginTop: theme.spacing(1),
    marginBottom: theme.spacing(1)
  }
}))

/**
 * Shows an informative overview about the selected entry.
 */
const OverviewView = React.memo(({entryId, ...moreProps}) => {
  const { raiseError } = useErrors()
  const [indexApiData, setIndexApiData] = useState(null)
  const [exists, setExists] = useState(true)
  const [archiveApiData, setArchiveApiData] = useState(null)
  const {api} = useApi()
  const index = useMemo(() => indexApiData?.response?.data, [indexApiData])
  const archive = useMemo(() => archiveApiData?.response?.data?.archive, [archiveApiData])
  const properties = useMemo(() => {
    return new Set(index?.results
      ? index.results.properties.available_properties
      : []
    )
  }, [index])

  useEffect(() => {
    async function fetchData() {
      const indexApiData = await api.get(`/entries/${entryId}`, null, {returnRequest: true})
      setIndexApiData(indexApiData)
      const archiveApiData = await api.post(
        `/entries/${entryId}/archive/query`,
        {
          required: {
            'resolve-inplace': false,
            results: {
              material: '*',
              method: '*',
              properties: {
                structures: '*',
                electronic: 'include-resolved',
                mechanical: 'include-resolved',
                spectroscopy: 'include-resolved',
                vibrational: 'include-resolved',
                // For geometry optimizations we require only the energies.
                // Trajectory, optimized structure, etc. are unnecessary.
                geometry_optimization: {
                  energies: 'include-resolved'
                }
              }
            }
          }
        },
        {returnRequest: true, jsonResponse: true}
      )
      setArchiveApiData(archiveApiData)
    }
    fetchData().catch(error => {
      if (error.name === 'DoesNotExist') {
        setExists(false)
      } else {
        raiseError(error)
      }
    })
  }, [api, raiseError, entryId, setIndexApiData, setExists, setArchiveApiData])

  const classes = useStyles()

  if (!exists) {
    return <Page>
      <Typography>
        This entry does not exist.
      </Typography>
    </Page>
  }

  if (!index) {
    return null
  }

  return <Page limitedWidth>
    <Grid container spacing={0} className={classes.root}>
      <Grid item xs={4} className={classes.leftColumn}>
        <MetadataSection title='Method'>
          <MethodMetadata data={index} />
        </MetadataSection>
        <Divider className={classes.divider} />
        <MetadataSection title='Author metadata'>
          <Quantity flex>
            <Quantity quantity='comment' data={index} />
            <Quantity quantity='references' data={index}/>
            <Quantity quantity='authors' data={index}/>
            <Quantity quantity="datasets" data={index}/>
          </Quantity>
        </MetadataSection>
        <Divider className={classes.divider}/>
        <MetadataSection>
          <Quantity column style={{maxWidth: 350}}>
            <Quantity quantity="mainfile" data={index}/>
            <Quantity quantity="entry_id" data={index}/>
            <Quantity quantity="results.material.material_id" data={index}/>
            <Quantity quantity="upload_id" data={index}/>
            <Quantity quantity="upload_create_time" data={index}/>
            <Quantity quantity="raw_id" data={index}/>
            <Quantity quantity="external_id" data={index}/>
            <Quantity quantity="last_processing_time" data={index}/>
            <Quantity quantity="last_processing_version" data={index}/>
          </Quantity>
        </MetadataSection>
        <SourceApiDialogButton label="API" maxWidth="lg" fullWidth buttonProps={{variant: 'contained', size: 'small'}}>
          {indexApiData && <SourceApiCall
            {...indexApiData}
            description="The basic metadata shown on this page is retrieved from the *entry metadata* API."
          />}
          <SourceDialogDivider />
          {archiveApiData && <SourceApiCall
            {...archiveApiData}
            description="The detailed property information is retrieved from the *entry archive* API. Only a specific parts of the archive are *required*."
          />}
        </SourceApiDialogButton>
      </Grid>

      <Grid item xs={8} className={classes.rightColumn}>
        <MaterialCard index={index} archive={archive} properties={properties}/>
        <ElectronicPropertiesCard index={index} archive={archive} properties={properties}/>
        <VibrationalPropertiesCard index={index} archive={archive} properties={properties}/>
        <MechanicalPropertiesCard index={index} archive={archive} properties={properties}/>
        <GeometryOptimizationCard index={index} archive={archive} properties={properties}/>
        <SpectroscopyCard index={index} archive={archive} properties={properties}/>
      </Grid>
    </Grid>
  </Page>
})

OverviewView.propTypes = {
  entryId: PropTypes.string.isRequired
}

OverviewView.whyDidYouRender = true

export default OverviewView
