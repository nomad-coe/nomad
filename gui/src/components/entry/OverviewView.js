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
import { Typography, makeStyles, Box, Grid, Link, Divider } from '@material-ui/core'
import { ApiDialog } from '../ApiDialogButton'
import { Actions } from '../Actions'
import Quantity from '../Quantity'
import { Link as RouterLink } from 'react-router-dom'
import { DOI } from '../dataset/DOI'
import { authorList } from '../../utils'
import searchQuantities from '../../searchQuantities'
import ElectronicPropertiesCard from '../entry/properties/ElectronicPropertiesCard'
import MaterialCard from '../entry/properties/MaterialCard'
import VibrationalPropertiesCard from '../entry/properties/VibrationalPropertiesCard'
import MechanicalPropertiesCard from '../entry/properties/MechanicalPropertiesCard'
import GeometryOptimizationCard from '../entry/properties/GeometryOptimizationCard'
import SpectroscopyCard from './properties/SpectroscopyCard'
import { MethodMetadata } from './EntryDetails'
import Page from '../Page'

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
const OverviewView = React.memo(function OverviewView({entryId, ...moreProps}) {
  const { raiseError } = useErrors()
  const [entry, setEntry] = useState(null)
  const [exists, setExists] = useState(true)
  const [showAPIDialog, setShowAPIDialog] = useState(false)
  const [archive, setArchive] = useState(null)
  const {api} = useApi()
  const properties = useMemo(() => {
    return new Set(entry?.results
      ? entry.results.properties.available_properties
      : []
    )
  }, [entry])

  useEffect(() => {
    api.entry(entryId).then(response => {
      const entry = response.data
      setEntry(entry)
      api.results(entry.entry_id)
        .then(setArchive)
        .catch(error => {
          if (error.name === 'DoesNotExist') {
          } else {
            raiseError(error)
          }
        })
    }).catch(error => {
      if (error.name === 'DoesNotExist') {
        setExists(false)
      } else {
        raiseError(error)
      }
    })
  }, [api, raiseError, entryId, setEntry, setExists, setArchive])

  const classes = useStyles()

  if (!exists) {
    return <Page>
      <Typography>
        This entry does not exist.
      </Typography>
    </Page>
  }

  if (!entry) {
    return null
  }

  return <Page limitedWidth>
    <Grid container spacing={0} className={classes.root}>
      <Grid item xs={4} className={classes.leftColumn}>
        <MetadataSection title='Method'>
          <MethodMetadata data={entry} />
        </MetadataSection>
        <Divider className={classes.divider} />
        <MetadataSection title='Author metadata'>
          <Quantity flex>
            <Quantity quantity='comment' placeholder='no comment' data={entry} />
            <Quantity quantity='references' placeholder='no references' data={entry}>
              {entry.references &&
              <div style={{display: 'inline-grid'}}>
                {entry.references.map(ref => <Typography key={ref} noWrap>
                  <Link href={ref}>{ref}</Link>
                </Typography>)}
              </div>}
            </Quantity>
            <Quantity quantity='authors' data={entry}>
              <Typography>
                {authorList(entry || [])}
              </Typography>
            </Quantity>
            <Quantity
              description={searchQuantities['datasets'] && searchQuantities['datasets'].description}
              label='datasets'
              placeholder='no datasets'
              data={entry}
            >
              {(entry.datasets && entry.datasets.length !== 0) &&
              <div>
                {entry.datasets.map(ds => (
                  <Typography key={ds.dataset_id}>
                    <Link component={RouterLink} to={`/dataset/id/${ds.dataset_id}`}>{ds.dataset_name}</Link>
                    {ds.doi ? <span>&nbsp;<DOI style={{display: 'inline'}} parentheses doi={ds.doi}/></span> : ''}
                  </Typography>))}
              </div>}
            </Quantity>
          </Quantity>
        </MetadataSection>
        <Divider className={classes.divider}/>
        <MetadataSection>
          <Quantity column style={{maxWidth: 350}}>
            <Quantity quantity="mainfile" noWrap ellipsisFront withClipboard data={entry}/>
            <Quantity quantity="entry_id" label='entry id' noWrap withClipboard data={entry}/>
            <Quantity quantity="results.material.material_id" label='material id' noWrap withClipboard data={entry}/>
            <Quantity quantity="upload_id" label='upload id' noWrap withClipboard data={entry}/>
            <Quantity quantity="upload_create_time" label='upload time' noWrap data={entry}>
              <Typography noWrap>
                {new Date(entry.upload_create_time).toLocaleString()}
              </Typography>
            </Quantity>
            <Quantity quantity="raw_id" label='raw id' noWrap hideIfUnavailable withClipboard data={entry}/>
            <Quantity quantity="external_id" label='external id' hideIfUnavailable noWrap withClipboard data={entry}/>
            <Quantity quantity="last_processing_time" label='last processing time' placeholder="not processed" noWrap data={entry}>
              <Typography noWrap>
                {new Date(entry.last_processing_time).toLocaleString()}
              </Typography>
            </Quantity>
            <Quantity description="Version used in the last processing" label='processing version' noWrap placeholder="not processed" data={entry}>
              <Typography noWrap>
                {entry.nomad_version}/{entry.nomad_commit}
              </Typography>
            </Quantity>
          </Quantity>
        </MetadataSection>
        <ApiDialog data={entry} open={showAPIDialog} onClose={() => { setShowAPIDialog(false) }}></ApiDialog>
        <Actions
          justifyContent='flex-end'
          variant='outlined'
          color='primary'
          size='medium'
          actions={[{
            tooltip: 'Show the API access code',
            onClick: (event) => { setShowAPIDialog(!showAPIDialog) },
            content: 'API'
          }]}
        >
        </Actions>
      </Grid>

      <Grid item xs={8} className={classes.rightColumn}>
        <MaterialCard entryMetadata={entry} archive={archive} properties={properties}/>
        <ElectronicPropertiesCard entryMetadata={entry} archive={archive} properties={properties}/>
        <VibrationalPropertiesCard entryMetadata={entry} archive={archive} properties={properties}/>
        <MechanicalPropertiesCard entryMetadata={entry} archive={archive} properties={properties}/>
        <GeometryOptimizationCard entryMetadata={entry} archive={archive} properties={properties}/>
        <SpectroscopyCard entryMetadata={entry} archive={archive} properties={properties}/>
      </Grid>
    </Grid>
  </Page>
})

OverviewView.propTypes = {
  entryId: PropTypes.string.isRequired
}

OverviewView.whyDidYouRender = true

export default OverviewView
