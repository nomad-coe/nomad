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
import 'regenerator-runtime/runtime'
import {render, screen, archives, wait, within} from '../../testSetup'
import { waitFor } from '@testing-library/dom'
import { getIndex } from '../../../tests/DFTBulk'
import { useApi } from '../api'
import { UploadsPage } from './UploadsPage'
import supertest from 'supertest'

const express = require('express')
const createMiddleware = require('@apidevtools/swagger-express-middleware')
let app = express()

// eslint-disable-next-line handle-callback-err
createMiddleware('openapi.json', app, (err, middleware) => {
  app.use(
    middleware.metadata(),
    middleware.parseRequest(),
    middleware.validateRequest(),
    middleware.mock()
  )

  supertest(app)
    .get('/info')
    .then(response => {
      if (response) {}
      // console.log(response.body)
    })
    // .end(function(err, res) {
    //   if (err) console.log(err)
    // })
})

jest.mock('../api')
const index = getIndex()

beforeAll(() => {
  // API mock init
  useApi.mockReturnValue({
    api: {
      post: () => wait({response: {data: {archive: archives.get(index.entry_id)}}}), // results
      get: (url) => {
        switch (url) {
          case '/uploads?page_size=10&page=1&order_by=upload_create_time&order=desc':
            return wait({response: { // response when page_size != 0
              query: {},
              data: [index],
              pagination: {
                'page_size': 10,
                'order_by': 'upload_create_time',
                'order': 'desc',
                'page': 1,
                'total': 1
              }
            }})
          case '/uploads?is_published=false&page_size=0&order_by=upload_create_time&order=desc':
            return wait(
              { // response when page_size=0
                query: {
                  'is_published': false
                },
                data: [],
                pagination: {
                  'page_size': 0,
                  'order_by': 'upload_create_time',
                  'order': 'desc',
                  'page': 1,
                  'total': 0
                }
              })
          default:
            return wait({response: { // response any other url
            }})
        }
      }}})
})

afterAll(() => {
  // API mock cleanup
  jest.unmock('../api')
})

test('correctly renders uploads page', async () => {
  render(<UploadsPage/>)

  // Wait to load the UploadsPage
  await waitFor(() => {
    expect(screen.queryByRole('new-upload-button')).toBeInTheDocument()
  })

  expect(screen.queryByRole('upload-commands')).toBeInTheDocument()
  expect(screen.queryByRole('error-maximum-number-of-unpublished')).toBeNull()

  // Test if the table header is rendered correctly
  expect(screen.queryByText('Your existing uploads')).toBeInTheDocument()
  expect(screen.queryByRole('table-pagination')).toBeInTheDocument()
  expect(screen.queryByRole('datatable-body')).toBeInTheDocument()

  let datatableBody = screen.getByRole('datatable-body')

  // Test if the times are printed correctly
  expect(within(datatableBody).queryByText('2021-03-17T13:47:32.899000')).toBeNull()
  expect(within(datatableBody).queryByText(new Date('2021-03-17T13:47:32.899000').toLocaleString())).toBeInTheDocument()

  // Test if the published icon is printed not unpublished icon
  expect(within(datatableBody).queryByTitle('published upload')).toBeInTheDocument()
  expect(within(datatableBody).queryByTitle('this upload is not yet published')).toBeNull()

  // TODO: test if the data is sorted by upload_create_time. It needs more test data in DFTBulk.js.
})
