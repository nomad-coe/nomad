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
import { React, useState } from 'react'

import { addSeconds, endOfDay, format, formatDistanceToNow, getUnixTime, isValid,
         subDays } from "date-fns"
import { Box, IconButton, Tooltip } from '@material-ui/core'
import ClipboardIcon from '@material-ui/icons/Assignment'
import { DatePicker } from '@material-ui/pickers'
import { CopyToClipboard } from 'react-copy-to-clipboard'
import { useKeycloak } from '@react-keycloak/web'

import { useApi } from './api'
import { appTokenMaxExpiresIn } from '../config'
import { useErrors } from './errors'

export default function AppTokenForm() {
  const {keycloak} = useKeycloak()
  const {api} = useApi()
  const errors = useErrors()

  const [error, setError] = useState()
  const [dateValue, setDateValue] = useState()
  const [token, setToken] = useState()

  let maxDate = subDays(addSeconds(new Date(), appTokenMaxExpiresIn), 1)
  if (!isValid(maxDate)) maxDate = undefined

  const handleDateChange = (newValue) => {
    setDateValue(newValue)
    if (error) setError(undefined)

    const expiresAt = endOfDay(newValue)
    const expiresInS = getUnixTime(expiresAt) - getUnixTime(new Date())
    api.get(`/auth/app_token?expires_in=${expiresInS}`)
      .then((response) => {
        setToken(response.app_token)
      })
      .catch((exception) => {
        setToken(undefined)

        if (exception.status !== 422) {
          errors.raiseError(exception)
        }

        const limitExpiresInS = exception.apiMessage[0].ctx.limit_value
        const limitExpiresAt = subDays(addSeconds(new Date(), limitExpiresInS), 1)
        setError(`Choose not later than ${format(limitExpiresAt, 'yyyy-MM-dd')} ` +
                 `(in ${formatDistanceToNow(limitExpiresAt)}).`)
      })
  }

  return (
    <Tooltip title={keycloak.authenticated ? "" : "Log in to request a token"}>
      <Box width='fit-content'>
        <DatePicker
          autoOk
          error={error}
          helperText={error}
          disabled={!keycloak.authenticated}
          disablePast
          maxDate={maxDate}
          size='small'
          variant='inline'
          inputVariant='filled'
          label='App token expires after'
          onChange={handleDateChange}
          format='yyyy-MM-dd'
          value={dateValue}
        />
        <CopyToClipboard text={token} >
          <Tooltip title="Copy token to clipboard">
            <IconButton disabled={!token}>
              <ClipboardIcon />
            </IconButton>
          </Tooltip>
        </CopyToClipboard>
      </Box>
    </Tooltip>
  )
}
