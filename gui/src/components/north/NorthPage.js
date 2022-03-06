/* eslint-disable quotes */
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
import React, { useState, useEffect, useCallback } from 'react'
import {
  Grid,
  Paper,
  Divider,
  Step,
  StepLabel,
  Stepper,
  makeStyles,
  Box
} from '@material-ui/core'
import Alert from '@material-ui/lab/Alert'
import Page from '../Page'
import { useApi, withLoginRequired } from '../api'
import { useErrors } from '../errors'
import NorthTool from './NorthTool'
import {
  addColumnDefaults,
  combinePagination,
  Datatable,
  DatatableLoadMorePagination,
  DatatableTable
} from '../datatable/Datatable'
import * as _tools from '../../northTools.json'

export const help = `
The NOMAD Remote Tools Hub (NORTH) provides access to tools which you can use to
work remotely on unpublished data stored within NOMAD.
`

// Datatable column setup
const columns = [
  {key: 'upload_id'},
  {key: 'upload_create_time'},
  {key: 'upload_name'}
]
addColumnDefaults(columns, {align: 'left'})

/**
 * Landing page for the NOMAD Remote Tools Hub.
 */

const useStyles = makeStyles(theme => ({
  stepper: {
    padding: 0,
    backgroundColor: 'inherit'
  },
  stepContent: {
    marginTop: theme.spacing(2),
    marginLeft: theme.spacing(3)
  }
}))
const NorthPage = React.memo(() => {
  const styles = useStyles()
  const tools = useTools()
  // const instances = useInstances()
  const {api} = useApi()
  const errors = useErrors()
  const [uploads, setUploads] = useState(null)
  const [selected, setSelected] = useState([])
  const [pagination, setPagination] = useState({
    page_size: 10,
    page: 1,
    order_by: 'upload_create_time',
    order: 'desc'
  })
  const canLaunch = selected.length > 0

  // Fetch the list of unpublished uploads.
  useEffect(() => {
    const {page_size, page, order_by, order} = pagination
    api.get(`/uploads?is_published=false&page_size=${page_size}&page=${page}&order_by=${order_by}&order=${order}`)
      .then(data => setUploads(data))
      .catch(errors.raiseError)
  }, [pagination, setUploads, errors, api])

  // Only one upload can be selected at this time
  const handleSelect = useCallback((callBack) => {
    const selected = callBack([])
    setSelected(old => {
      if (old[0]?.upload_id === selected[0]?.upload_id) {
        return []
      }
      return selected
    })
  }, [])

  return (uploads)
    ? <Page limitedWidth>
      <Grid container spacing={1}>
        <Grid item xs={8}>
          <Stepper className={styles.stepper} orientation="horizontal">
            <Step active>
              <StepLabel icon={1}>Select an upload</StepLabel>
            </Step>
          </Stepper>
          <Paper className={styles.stepContent}>
            {uploads.data.length > 0
              ? <Datatable
                multiple={false}
                selected={selected}
                onSelectedChanged={handleSelect}
                columns={columns}
                selectedColumns={columns.map(column => column.key)}
                data={uploads.data || []}
                pagination={combinePagination(pagination, uploads.pagination)}
                onPaginationChanged={setPagination}
              >
                <DatatableTable>
                  <DatatableLoadMorePagination />
                </DatatableTable>
              </Datatable>
              : <Alert severity="info">No uploads available. Notice that you can only work on unpublished uploads.</Alert>
            }
          </Paper>
        </Grid>
        <Grid item xs={4}>
          <Stepper className={styles.stepper} orientation="horizontal">
            <Step active={canLaunch}>
              <StepLabel icon={2}>Launch a tool</StepLabel>
            </Step>
          </Stepper>
          <Paper className={styles.stepContent}>
            {Object.keys(tools).map(key => ({name: key, title: key, ...tools[key]})).map((tool, index) => (
              <React.Fragment key={tool.name}>
                <Box padding={1}>
                  <NorthTool
                    {...tool}
                    disabled={!canLaunch}
                    uploadId={selected[0] && selected[0].upload_id}
                  />
                </Box>
                {index !== Object.keys(tools).length - 1 && (
                  <Box marginY={1}><Divider/></Box>
                )}
              </React.Fragment>
            ))}
          </Paper>
        </Grid>
      </Grid>
    </Page>
    : null
})

export default withLoginRequired(NorthPage)

/**
 * Hook for loading the list of running instances from the NORTH API.
*/
// TODO use Jupyterhub for this
export function useInstances() {
  const [instances, setInstances] = useState()
  const { api } = useApi()

  useEffect(() => {
    api.northInstances()
      .then(data => {
        const instances = {}
        data.forEach(instance => {
          instances[instance.name] = {...instance}
        })
        setInstances(instances)
      })
  }, [api])

  return instances
}

/**
 * Hook for loading the list of available tools from the NORTH API.
*/
export function useTools() {
  return _tools.default
}
