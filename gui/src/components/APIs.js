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
import { apiBase, appBase } from '../config'
import Markdown from './Markdown'
import AppTokenForm from './AppTokenForm'

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
      and functions programatically. For all APIs, we offer dashboards that let you use
      each API interactively, right in your browser.

      ## NOMAD's API

      - [API dashboard](${apiBase}/v1/extensions/docs)
      - [API documentation](${apiBase}/v1/extensions/redoc)

      We started to implement a more consise and easier to use API for access NOMAD
      data. This will step-by-step reimplement all functions of NOMAD's old main API.
      At some point, it will replace it entirely. For new users, we recommend to start
      using this API. API Dashboard and documentation contain a tutorial on how to get started.

      There is a [tutorial on how to use the API with plain Python](${appBase}/docs/api_tutorial.html).
      Another [tutorial covers how to install and use NOMAD's Python client library](${appBase}/docs/archive_tutorial.html).
      The [NOMAD Analytics Toolkit](https://nomad-lab.eu/AIToolkit) allows to use
      this without installation and directly on NOMAD servers.

      ###  App token

      Next to the usual access token based on OpenID Connect, we provide an
      [app token](${appBase}/docs/apis/api.html#app-token) with custom expiration date.
      You may request one via the [API](${apiBase}/v1/extensions/docs) or the following form:
    `}</Markdown>

    <AppTokenForm />

    <Markdown>{`

      ### Old API

      You can still use NOMAD's old REST API. The data it provides might miss the most
      recent contributions:

      - [v0 API dashboard](https://nomad-lab.eu/prod/rae/api/)

      ## OPTIMADE

      - [OPTIMADE API overview page](${appBase}/optimade/)
      - [OPTIMADE API dashboard](${appBase}/optimade/v1/extensions/docs)
      - [OPTIMADE API documentation](${appBase}/optimade/v1/extensions/redoc)

      [OPTIMADE](https://www.optimade.org/) is an
      open API standard for materials science databases. This API can be used to search
      and access NOMAD metadata in a standardized way that can also be applied to many
      [other materials science databses](https://providers.optimade.org/).

      ## DCAT

      - [DCAT API dashboard](${appBase}/dcat/)

      [DCAT](https://www.w3.org/TR/vocab-dcat-2/) is a RDF vocabulary designed to facilitate
      interoperability between data catalogs published on the Web. This API allows you
      access to NOMAD via RDF documents following DCAT. You can access NOMAD entries as
      DCAT Datasets or all NOMAD entries as a DCAT Catalog.

      ## Resources

      - [Resources API dashboard](${appBase}/resources/extensions/docs)

      The resources API provides links from NOMAD entries to related external resources.
      These include the [Aflow Encyclopedia of Crystallographic Prototypes](https://www.aflowlib.org/prototype-encyclopedia/),
      [Springer Materials Database of Inorganic Solid Phases](https://materials.springer.com)
      and [OPTIMADE providers](https://providers.optimade.org/).
    `}</Markdown>
  </div>
}
