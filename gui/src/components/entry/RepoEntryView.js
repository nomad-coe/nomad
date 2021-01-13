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
import React, { useContext, useState, useEffect } from 'react'
import PropTypes from 'prop-types'
import { Divider, Card, CardContent, Grid, CardHeader, Typography, Link, makeStyles } from '@material-ui/core'
import { apiContext } from '../api'
import ApiDialogButton from '../ApiDialogButton'
// import Structure from '../visualization/Structure'
import Quantity from '../Quantity'
import { Link as RouterLink } from 'react-router-dom'
import { DOI } from '../search/DatasetList'
import { domains } from '../domains'
import { EntryPageContent } from './EntryPage'
import { errorContext } from '../errors'
import { authorList } from '../../utils'

const useStyles = makeStyles(theme => ({
  root: {
    marginTop: theme.spacing(2)
  },
  error: {
    marginTop: theme.spacing(2)
  },
  cardContent: {
    paddingTop: 0
  },
  entryCards: {
    marginTop: theme.spacing(2)
  },
  structureViewer: {
    height: '25rem',
    padding: '0'
  }
}))

export default function RepoEntryView({uploadId, calcId}) {
  const classes = useStyles()
  const {api} = useContext(apiContext)
  const {raiseError} = useContext(errorContext)
  const [state, setState] = useState({calcData: null, doesNotExist: false})

  useEffect(() => {
    setState({calcData: null, doesNotExist: false})
  }, [setState, uploadId, calcId])

  useEffect(() => {
    api.repo(uploadId, calcId).then(data => {
      setState({calcData: data, doesNotExist: false})
    }).catch(error => {
      if (error.name === 'DoesNotExist') {
        setState({calcData: null, doesNotExist: true})
      } else {
        setState({calcData: null, doesNotExist: false})
        raiseError(error)
      }
    })
  }, [api, raiseError, uploadId, calcId, setState])

  const calcData = state.calcData || {uploadId: uploadId, calcId: calcId}
  const loading = !state.calcData
  const quantityProps = {data: calcData, loading: loading}

  const domain = calcData.domain && domains[calcData.domain]

  let entryHeader = 'Entry metadata'
  if (domain) {
    entryHeader = domain.entryTitle(calcData)
  }

  if (state.doesNotExist) {
    return <EntryPageContent className={classes.root} fixed>
      <Typography className={classes.error}>
        This entry does not exist.
      </Typography>
    </EntryPageContent>
  }

  return (
    <EntryPageContent className={classes.root} fixed>
      <Grid container spacing={2}>
        <Grid item xs={7}>
          <Card>
            <CardHeader
              title={entryHeader}
              action={<ApiDialogButton title="Repository JSON" data={calcData} />}
            />
            <CardContent classes={{root: classes.cardContent}}>
              {domain && <domain.EntryOverview data={calcData} loading={loading} />}
            </CardContent>
            <Divider />
            <CardContent>
              <Quantity column>
                <Quantity quantity='comment' placeholder='no comment' {...quantityProps} />
                <Quantity quantity='references' placeholder='no references' {...quantityProps}>
                  {calcData.references &&
                    <div style={{display: 'inline-grid'}}>
                      {calcData.references.map(ref => <Typography key={ref} noWrap>
                        <Link href={ref}>{ref}</Link>
                      </Typography>)}
                    </div>}
                </Quantity>
                <Quantity quantity='authors' {...quantityProps}>
                  <Typography>
                    {authorList(loading ? null : calcData, true)}
                  </Typography>
                </Quantity>
                <Quantity quantity='datasets' placeholder='no datasets' {...quantityProps}>
                  {calcData.datasets &&
                    <div>
                      {calcData.datasets.map(ds => (
                        <Typography key={ds.dataset_id}>
                          <Link component={RouterLink} to={`/dataset/id/${ds.dataset_id}`}>{ds.name}</Link>
                          {ds.doi ? <span>&nbsp; (<DOI doi={ds.doi}/>)</span> : ''}
                        </Typography>))}
                    </div>}
                </Quantity>
              </Quantity>
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={5}>
          <Card>
            <CardHeader title="Ids / processing" />
            <CardContent classes={{root: classes.cardContent}}>
              <Quantity column style={{maxWidth: 350}}>
                <Quantity quantity="calc_id" label={`${domain ? domain.entryLabel : 'entry'} id`} noWrap withClipboard {...quantityProps} />
                <Quantity quantity="encyclopedia.material.material_id" label='material id' loading={loading} noWrap {...quantityProps} withClipboard />
                <Quantity quantity="mainfile" loading={loading} noWrap ellipsisFront {...quantityProps} withClipboard />
                <Quantity quantity="upload_id" label='upload id' {...quantityProps} noWrap withClipboard />
                <Quantity quantity="upload_time" label='upload time' noWrap {...quantityProps} >
                  <Typography noWrap>
                    {new Date(calcData.upload_time).toLocaleString()}
                  </Typography>
                </Quantity>
                <Quantity quantity="raw_id" label='raw id' loading={loading} noWrap hideIfUnavailable {...quantityProps} withClipboard />
                <Quantity quantity="external_id" label='external id' loading={loading} hideIfUnavailable noWrap {...quantityProps} withClipboard />
                <Quantity quantity="last_processing" label='last processing' loading={loading} placeholder="not processed" noWrap {...quantityProps}>
                  <Typography noWrap>
                    {new Date(calcData.last_processing).toLocaleString()}
                  </Typography>
                </Quantity>
                <Quantity quantity="last_processing" label='processing version' loading={loading} noWrap placeholder="not processed" {...quantityProps}>
                  <Typography noWrap>
                    {calcData.nomad_version}/{calcData.nomad_commit}
                  </Typography>
                </Quantity>
              </Quantity>
            </CardContent>
          </Card>
        </Grid>
      </Grid>

      {domain && <domain.EntryCards data={calcData} calcId={calcId} uploadId={uploadId} classes={{root: classes.entryCards}} />}
    </EntryPageContent>
  )
}

RepoEntryView.propTypes = {
  uploadId: PropTypes.string.isRequired,
  calcId: PropTypes.string.isRequired
}
