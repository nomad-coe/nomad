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
import { makeStyles } from '@material-ui/core'
import React from 'react'
import { apiBase, appBase, optimadeBase } from '../config'
import Markdown from './Markdown'

const useStyles = makeStyles(theme => ({
  root: {
    padding: theme.spacing(3),
    maxWidth: 1024,
    margin: 'auto',
    width: '100%'
  }
}))

export default function About() {
  const classes = useStyles()

  return <div className={classes.root}>
    <Markdown>{`
      # APIs

      NOMAD's Application Programming Interface (API) allows you to access NOMAD data
      and functions programatically.

      ## NOMAD's main API

      - [API dashboard](${apiBase}/)

      This is NOMAD main REST API. This API the main interface to NOMAD and it also used
      by this web-page to provide all functions. Therefore, everything you do here, can
      also be done by using this API.

      There is a [tutorial on how to use the API with plain Python](${appBase}/docs/api_tutorial.html).
      Another [tutorial covers how to install and use NOMAD's Python client library](${appBase}/docs/archive_tutorial.html).
      The [NOMAD Analytics Toolkit](https://nomad-lab.eu/AIToolkit) allows to use
      this without installation and directly on NOMAD servers.

      ## OPTIMADE

      - [OPTIMADE API dashboard](${optimadeBase}/)

      [OPTIMADE](https://www.optimade.org/) is an
      open API standard for materials science databases. This API can be used to search
      and access NOMAD metadata in a standardized way that can also be applied to many
      [other materials science databses](https://providers.optimade.org/).
    `}</Markdown>
  </div>
}
