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
import React, { useEffect, useMemo, useState } from 'react'
import PropTypes from 'prop-types'
import { IconButton, Paper, Tabs, Tab, Box } from '@material-ui/core'
import TooltipButton from '../utils/TooltipButton'
import EditIcon from '@material-ui/icons/Edit'
import NavigateNextIcon from '@material-ui/icons/NavigateNext'
import { Datatable, DatatableToolbar, DatatableToolbarActions, DatatableTable,
  DatatablePagePagination, DatatableScrollPagination, DatatableLoadMorePagination } from './Datatable'
import Page from '../Page'
import { sleep } from '../../utils'

const createExampleData = (length) => Array.from({length: length}, (_, i) => ({
  id: (i + 1).toString(), value: i + 1
}))

const columns = [
  {
    key: 'id',
    sortable: true,
    align: 'right',
    render: data => data.id
  },
  {
    key: 'value',
    sortable: true,
    render: data => `This is value for id ${data.value}.`
  }
]

function DetailsExample({data}) {
  return <div>Details: {data.id}</div>
}
DetailsExample.propTypes = {
  data: PropTypes.object.isRequired
}

function ActionsExample({data}) {
  return <TooltipButton
    title="A row action." component={IconButton}
    onClick={() => console.log(`Perform row action for ${data.id}.`)}
  >
    <NavigateNextIcon/>
  </TooltipButton>
}
ActionsExample.propTypes = {
  data: PropTypes.object.isRequired
}

function ScrollDatatableExample() {
  const exampleData = useMemo(() => createExampleData(1024), [])
  const [pagination, setPagination] = useState({
    total: 1024, page_size: 100
  })
  const [data, setData] = useState([])
  const [selected, setSelected] = useState(new Set())

  useEffect(() => {
    sleep(100).then(() => {
      const nextPageAfterValueInt = parseInt(pagination.page_after_value || 0) + pagination.page_size
      setData(exampleData.slice(0, nextPageAfterValueInt))
      setPagination(pagination => ({
        ...pagination,
        next_page_after_value: nextPageAfterValueInt > exampleData.length ? null : String(nextPageAfterValueInt)
      }))
    })
  }, [pagination.page_after_value, pagination.page_size, exampleData])

  return <Paper style={{height: '100%'}}>
    <Datatable
      columns={columns} data={data}
      pagination={pagination} onPaginationChanged={setPagination}
      selected={selected}
      getId={option => option.id}
      onSelectedChanged={setSelected}
    >
      <DatatableToolbar title="Example table." />
      <DatatableTable
        details={DetailsExample}
        actions={ActionsExample}
      >
        <DatatableScrollPagination loadingMessage="loading more example items ..."/>
      </DatatableTable>
    </Datatable>
  </Paper>
}

function LoadMoreDatatableExample() {
  const exampleData = useMemo(() => createExampleData(1024), [])
  const [pagination, setPagination] = useState({
    total: 1024, page_size: 10
  })
  const [data, setData] = useState([])
  const [selected, setSelected] = useState(new Set())

  useEffect(() => {
    sleep(100).then(() => {
      const nextPageAfterValueInt = parseInt(pagination.page_after_value || 0) + pagination.page_size
      setData(exampleData.slice(0, nextPageAfterValueInt))
      setPagination(pagination => ({
        ...pagination,
        next_page_after_value: nextPageAfterValueInt > exampleData.length ? null : String(nextPageAfterValueInt)
      }))
    })
  }, [pagination.page_after_value, pagination.page_size, exampleData])

  return <Paper>
    <Datatable
      columns={columns} data={data}
      pagination={pagination} onPaginationChanged={setPagination}
      selected={selected}
      getId={option => option.id}
      onSelectedChanged={setSelected}
    >
      <DatatableToolbar title="Example table." />
      <DatatableTable
        details={DetailsExample}
        actions={ActionsExample}
      >
        <DatatableLoadMorePagination color="primary">
          load more
        </DatatableLoadMorePagination>
      </DatatableTable>
    </Datatable>
  </Paper>
}

export function DatatableExamples() {
  const exampleData = useMemo(() => createExampleData(23), [])

  const [data, setData] = useState([])
  const [pagination, setPagination] = useState({
    page_size: 10, page: 1, order_by: 'id', order: 'asc', total: exampleData.length
  })
  const [selected, setSelected] = useState(new Set())

  useEffect(() => {
    const {page_size, page, total, order, order_by} = pagination
    const compare = (a, b) => {
      const value = data => parseInt(String(data[order_by]))
      return (value(a) - value(b)) * (order === 'asc' ? 1 : -1)
    }
    const data = [...exampleData].sort(compare)
    const start = page_size * (page - 1)
    const end = page_size + start
    const length = (end > total) ? total - start : page_size
    setData(data.slice(start, start + length))
  }, [pagination, exampleData])

  const [tab, setTab] = useState(3)

  return <Box height="100%" display="flex" flexDirection="column">
    <Tabs value={tab} onChange={(event, tab) => setTab(tab)}>
      <Tab label="details paginated" />
      <Tab label="paginated" />
      <Tab label="plain" />
      <Tab label="scroll" />
      <Tab label="load more" />
    </Tabs>

    <Box flex={1}>
      <Page limitedWidth hidden={tab !== 0}>
        <Paper>
          <Datatable
            columns={columns}
            data={data}
            pagination={pagination}
            onPaginationChanged={setPagination}
            selected={selected}
            getId={option => option.id}
            onSelectedChanged={setSelected}
          >
            <DatatableToolbar title="Example table.">
              <DatatableToolbarActions selection>
                <IconButton onClick={() => console.log('Perform action on selection.')}>
                  <EditIcon/>
                </IconButton>
              </DatatableToolbarActions>
              <DatatableToolbarActions>
                <IconButton onClick={() => console.log('Perform general action.')}>
                  <EditIcon/>
                </IconButton>
              </DatatableToolbarActions>
            </DatatableToolbar>
            <DatatableTable details={DetailsExample} actions={ActionsExample}>
              <DatatablePagePagination pageSizeValues={[5, 10, 50, 100]}/>
            </DatatableTable>
          </Datatable>
        </Paper>
      </Page>

      <Page limitedWidth hidden={tab !== 1}>
        <Paper>
          <Datatable
            columns={columns}
            data={data}
            pagination={pagination}
            onPaginationChanged={setPagination}
            selected={selected}
            getId={option => option.id}
            onSelectedChanged={setSelected}
          >
            <DatatableToolbar title="Example table.">
              <DatatableToolbarActions selection>
                <IconButton onClick={() => console.log('Perform action on selection.')}>
                  <EditIcon/>
                </IconButton>
              </DatatableToolbarActions>
              <DatatableToolbarActions>
                <IconButton onClick={() => console.log('Perform general action.')}>
                  <EditIcon/>
                </IconButton>
              </DatatableToolbarActions>
            </DatatableToolbar>
            <DatatableTable>
              <DatatablePagePagination pageSizeValues={[5, 10, 50, 100]}/>
            </DatatableTable>
          </Datatable>
        </Paper>
      </Page>
      <Page limitedWidth hidden={tab !== 2}>
        <Paper>
          <Datatable columns={columns} data={exampleData} />
        </Paper>
      </Page>
      <Page limitedWidth limitedHeight hidden={tab !== 3}>
        <ScrollDatatableExample/>
      </Page>
      <Page limitedWidth limitedHeight hidden={tab !== 4}>
        <LoadMoreDatatableExample/>
      </Page>
    </Box>
  </Box>
}
