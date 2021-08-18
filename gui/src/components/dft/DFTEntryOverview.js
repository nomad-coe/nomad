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
import { Box, CardContent, Grid, Typography, Link, makeStyles, Divider } from '@material-ui/core'
import { useApi } from '../apiV1'
import { ApiDialog } from '../ApiDialogButton'
import Actions from '../Actions'
import Quantity from '../Quantity'
import { Link as RouterLink } from 'react-router-dom'
import { DOI } from '../search/DatasetList'
import { useErrors } from '../errors'
import { authorList } from '../../utils'
import searchQuantities from '../../searchQuantities'
import ElectronicPropertiesCard from '../entry/properties/ElectronicPropertiesCard'
import MaterialCard from '../entry/properties/MaterialCard'
import VibrationalPropertiesCard from '../entry/properties/VibrationalPropertiesCard'
import GeometryOptimizationCard from '../entry/properties/GeometryOptimizationCard'

const useHeaderStyles = makeStyles(theme => ({
  root: {
    marginBottom: theme.spacing(1),
    display: 'flex',
    flexDirection: 'row',
    alignContent: 'flex-start'
  },
  title: {
    fontSize: '1.25rem'
  }
}))

function Header({title, actions}) {
  const styles = useHeaderStyles()
  return <Box className={styles.root}>
    <Box flexGrow={1}>
      <Typography className={styles.title}>
        {title}
      </Typography>
    </Box>
    {actions}
  </Box>
}

Header.propTypes = {
  title: PropTypes.string,
  actions: PropTypes.any
}

const useSidebarCardStyles = makeStyles(theme => ({
  content: {
    padding: 0,
    paddingBottom: theme.spacing(2)
  },
  title: {
    fontSize: '1.25rem',
    marginBottom: theme.spacing(1)
  }
}))

function SidebarCard({title, actions, children}) {
  const styles = useSidebarCardStyles()
  return <CardContent className={styles.content}>
    {(title || actions) && <Header title={title} actions={actions}></Header>}
    {children}
  </CardContent>
}

SidebarCard.propTypes = {
  title: PropTypes.string,
  actions: PropTypes.any,
  children: PropTypes.any
}

const useStyles = makeStyles(theme => ({
  root: {
    marginBottom: theme.spacing(4)
  },
  error: {
    marginTop: theme.spacing(2)
  },
  leftSidebar: {
    maxWidth: '32%',
    flexBasis: '32%',
    flexGrow: 0,
    paddingRight: theme.spacing(3)
  },
  rightSidebar: {
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
const DFTEntryOverview = ({data}) => {
  const api = useApi()
  const {raiseError} = useErrors()
  const [showAPIDialog, setShowAPIDialog] = useState(false)
  const [archive, setArchive] = useState(null)

  const styles = useStyles()

  // Load the archive data from the server. If section_results is present, only
  // part of the archive will be loaded. Otherwise a full archive download is
  // performed.
  useEffect(() => {
    api.results(data.entry_id)
      .then(setArchive)
      .catch(error => {
        if (error.name === 'DoesNotExist') {
        } else {
          raiseError(error)
        }
      })
  }, [data.entry_id, api, raiseError, setArchive])

  const methodQuantities = useMemo(() => {
    const methodQuantities = []
    const addMethodQuantities = (obj, parentKey) => {
      const children = {}
      Object.keys(obj).forEach(key => {
        const value = obj[key]
        if (Array.isArray(value) || typeof value === 'string') {
          if (value.length > 0) {
            methodQuantities.push({
              quantity: `${parentKey}.${key}`,
              label: key.replaceAll('_', ' ')
            })
          }
        } else if (value instanceof Object) {
          children[key] = value
        }
      })
      Object.keys(children).forEach(key => addMethodQuantities(children[key], `${parentKey}.${key}`))
    }
    addMethodQuantities(data.results.method, 'results.method')
    return methodQuantities
  }, [data])

  return (
    <Grid container spacing={0} className={styles.root}>
      {/* Left column */}
      <Grid item xs={4} className={styles.leftSidebar}>
        <SidebarCard title='Method'>
          <Quantity flex>
            {methodQuantities.map(({...quantityProps}) => (
              <Quantity
                key={quantityProps.quantity}
                {...quantityProps}
                noWrap
                data={data}
                hideIfUnavailable
              />
            ))}
          </Quantity>
        </SidebarCard>
        <Divider className={styles.divider} />
        <SidebarCard title='Author metadata'>
          <Quantity flex>
            <Quantity quantity='comment' placeholder='no comment' data={data} />
            <Quantity quantity='references' placeholder='no references' data={data}>
              {data.references &&
              <div style={{display: 'inline-grid'}}>
                {data.references.map(ref => <Typography key={ref} noWrap>
                  <Link href={ref}>{ref}</Link>
                </Typography>)}
              </div>}
            </Quantity>
            <Quantity quantity='authors' data={data}>
              <Typography>
                {authorList(data || [])}
              </Typography>
            </Quantity>
            <Quantity
              description={searchQuantities['datasets'] && searchQuantities['datasets'].description}
              label='datasets'
              placeholder='no datasets'
              data={data}
            >
              {(data.datasets && data.datasets.length !== 0) &&
              <div>
                {data.datasets.map(ds => (
                  <Typography key={ds.dataset_id}>
                    <Link component={RouterLink} to={`/dataset/id/${ds.dataset_id}`}>{ds.name}</Link>
                    {ds.doi ? <span>&nbsp;<DOI style={{display: 'inline'}} parentheses doi={ds.doi}/></span> : ''}
                  </Typography>))}
              </div>}
            </Quantity>
          </Quantity>
        </SidebarCard>
        <Divider className={styles.divider}/>
        <SidebarCard>
          <Quantity column style={{maxWidth: 350}}>
            <Quantity quantity="mainfile" noWrap ellipsisFront withClipboard data={data}/>
            <Quantity quantity="entry_id" label='entry id' noWrap withClipboard data={data}/>
            <Quantity quantity="results.material.material_id" label='material id' noWrap withClipboard data={data}/>
            <Quantity quantity="upload_id" label='upload id' noWrap withClipboard data={data}/>
            <Quantity quantity="upload_time" label='upload time' noWrap data={data}>
              <Typography noWrap>
                {new Date(data.upload_time).toLocaleString()}
              </Typography>
            </Quantity>
            <Quantity quantity="raw_id" label='raw id' noWrap hideIfUnavailable withClipboard data={data}/>
            <Quantity quantity="external_id" label='external id' hideIfUnavailable noWrap withClipboard data={data}/>
            <Quantity quantity="last_processing" label='last processing' placeholder="not processed" noWrap data={data}>
              <Typography noWrap>
                {new Date(data.last_processing).toLocaleString()}
              </Typography>
            </Quantity>
            <Quantity description="Version used in the last processing" label='processing version' noWrap placeholder="not processed" data={data}>
              <Typography noWrap>
                {data.nomad_version}/{data.nomad_commit}
              </Typography>
            </Quantity>
          </Quantity>
        </SidebarCard>
        <ApiDialog data={data} open={showAPIDialog} onClose={() => { setShowAPIDialog(false) }}></ApiDialog>
        <Actions justifyContent='flex-end'>
          <Action
            tooltip="Show the API access code"
            onClick={(event) => { setShowAPIDialog(!showAPIDialog) }}
            variant="outlined"
            color="primary"
            size="medium"
          >
            API
          </Action>
        </Actions>
      </Grid>

      {/* Right column */}
      <Grid item xs={8} className={styles.rightSidebar}>
        <MaterialCard entryMetadata={data} archive={archive} />
        <ElectronicPropertiesCard entryMetadata={data} archive={archive} />
        <VibrationalPropertiesCard entryMetadata={data} archive={archive} />
        <GeometryOptimizationCard entryMetadata={data} archive={archive} />
      </Grid>
    </Grid>
  )
}

DFTEntryOverview.propTypes = {
  data: PropTypes.object.isRequired
}

export default DFTEntryOverview
