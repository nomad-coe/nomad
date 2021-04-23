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
import Quantity from '../Quantity'
import { Divider, Grid, Card, CardHeader, CardContent, Typography, Link, makeStyles } from '@material-ui/core'
import { domainData } from '../domainData'
import ApiDialogButton from '../ApiDialogButton'
import { Link as RouterLink } from 'react-router-dom'
import { DOI } from '../search/DatasetList'
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

export default function DefaultEntryOverview({data, children}) {
  const classes = useStyles()
  const domain = data.domain && domainData[data.domain]
  let entryHeader = 'Entry metadata'
  if (domain) {
    entryHeader = domain.entryTitle(data)
  }

  return (
    <Grid container spacing={2}>
      <Grid item xs={7}>
        <Card>
          <CardHeader
            title={entryHeader}
            action={<ApiDialogButton title="Repository JSON" data={data} />}
          />
          <CardContent classes={{root: classes.cardContent}}>
            {children}
          </CardContent>
          <Divider />
          <CardContent>
            <Quantity column>
              <Quantity quantity='comment' placeholder='no comment' data={data}/>
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
                  {authorList(data, true)}
                </Typography>
              </Quantity>
              <Quantity quantity='datasets' data={data}>
                {data.datasets && data.datasets.length > 0
                  ? <div>
                    {data.datasets.map(ds => (
                      <Typography key={ds.dataset_id}>
                        <Link component={RouterLink} to={`/dataset/id/${ds.dataset_id}`}>{ds.name}</Link>
                        {ds.doi ? <span>&nbsp; (<DOI doi={ds.doi} />)</span> : ''}
                      </Typography>))}
                  </div>
                  : <Typography><i>not in any dataset</i></Typography>}
              </Quantity>
              <Quantity quantity='license' placeholder='unspecified' data={data}/>
            </Quantity>
          </CardContent>
        </Card>
      </Grid>

      <Grid item xs={5}>
        <Card>
          <CardHeader title="Ids / processing" />
          <CardContent classes={{root: classes.cardContent}}>
            <Quantity column style={{maxWidth: 350}}>
              <Quantity quantity="calc_id" label={`${domain ? domain.entryLabel : 'entry'} id`} noWrap withClipboard data={data} />
              <Quantity quantity="encyclopedia.material.material_id" label='material id' noWrap data={data} withClipboard />
              <Quantity quantity="mainfile" noWrap ellipsisFront data={data} withClipboard />
              <Quantity quantity="upload_id" label='upload id' data={data} noWrap withClipboard />
              <Quantity quantity="upload_time" label='upload time' noWrap data={data}>
                <Typography noWrap>
                  {new Date(data.upload_time).toLocaleString()}
                </Typography>
              </Quantity>
              <Quantity quantity="raw_id" label='raw id' noWrap hideIfUnavailable data={data} withClipboard />
              <Quantity quantity="external_id" label='external id' hideIfUnavailable noWrap data={data} withClipboard />
              <Quantity quantity="last_processing" label='last processing' placeholder="not processed" noWrap data={data}>
                <Typography noWrap>
                  {new Date(data.last_processing).toLocaleString()}
                </Typography>
              </Quantity>
              <Quantity quantity="last_processing" label='processing version' noWrap placeholder="not processed" dat={data}>
                <Typography noWrap>
                  {data.nomad_version}/{data.nomad_commit}
                </Typography>
              </Quantity>
            </Quantity>
          </CardContent>
        </Card>
      </Grid>
    </Grid>
  )
}

DefaultEntryOverview.propTypes = {
  data: PropTypes.object.isRequired,
  children: PropTypes.object
}
