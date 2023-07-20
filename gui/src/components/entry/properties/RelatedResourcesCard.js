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
import { Link, IconButton, Typography, makeStyles, Box } from '@material-ui/core'
import TooltipButton from '../../utils/TooltipButton'
import NavigateNextIcon from '@material-ui/icons/NavigateNext'
import { useApi } from '../../api'
import { useErrors } from '../../errors'
import { Datatable, DatatableTable, DatatablePagePagination } from '../../datatable/Datatable'
import { formatTimestamp } from '../../../utils'
import Quantity from '../../Quantity'
import { PropertyCard } from './PropertyCard'
import Ellipsis from '../../visualization/Ellipsis'

const useResourceDetailsStyles = makeStyles(theme => ({
  resourceDetails: {
    padding: theme.spacing(1)
  },
  resourceDetailsContents: {
    display: 'flex',
    width: '100%',
    margin: '0'
  },
  resourceDetailsRow: {
    paddingRight: theme.spacing(3)
  },
  resourceDetailsActions: {
    display: 'flex',
    flexBasis: 'auto',
    flexGrow: 0,
    flexShrink: 0,
    justifyContent: 'flex-end',
    marginBottom: theme.spacing(1),
    marginTop: theme.spacing(2)
  },
  resourceURL: {
    overflow: 'hidden',
    textOverflow: 'ellipsis',
    width: '11rem'
  }
}))

const ResourceDetails = React.memo(({data}) => {
  const classes = useResourceDetailsStyles()
  const downloadTime = {download_time: formatTimestamp(data?.download_time)}

  return <div className={classes.resourceDetails}>
    <div className={classes.resourceDetailsContents}>
      <div className={classes.resourceDetailsRow}>
        <Quantity quantity='download_time' label='Last time accessed' data={downloadTime}></Quantity>
        <Quantity quantity='available_data' label='Available data' data={data}></Quantity>
      </div>
      <div className={classes.resourceDetailsRow}>
        <Quantity quantity='database_version' label='Database version' data={data}></Quantity>
      </div>
      <div className={classes.resourceDetailsRow}>
        <Quantity quantity='comment' label='Comment' data={data}></Quantity>
      </div>
    </div>
  </div>
})

ResourceDetails.propTypes = {
  data: PropTypes.object.isRequired
}

const ResourceActions = React.memo(({data}) => {
  return <TooltipButton
    title="Go to resource page." component={IconButton}
    onClick={() => { window.location = data.url }}
  >
    <NavigateNextIcon/>
  </TooltipButton>
})

ResourceActions.propTypes = {
  data: PropTypes.object.isRequired
}

const RelatedResourcesCard = React.memo(({index, archive}) => {
  const {resourcesApi} = useApi()
  const {raiseError} = useErrors()
  const [externalResources, setExternalResources] = useState()
  const [pagination, setPagination] = useState({
    page_size: 10, page: 1, order_by: 'id', order: 'asc', total: 0
  })

  useEffect(() => {
    const material = archive?.results?.material
    const spg = material?.symmetry?.space_group_number
    const structures = archive?.results?.properties?.structures
    const wyckoffSets = structures?.structure_conventional?.wyckoff_sets
    const wyckoff = wyckoffSets ? wyckoffSets.map((set) => set?.wyckoff_letter || '') : null
    const species = structures?.structure_original?.species_at_sites || null
    const nsites = species ? species.length : null
    const formula = material?.chemical_formula_reduced
    const params = []
    if (spg) {
      params.push(`space_group_number=${spg}`)
    }
    if (wyckoff) {
      params.push(`wyckoff_letters=${wyckoff.join('&wyckoff_letters=')}`)
    }
    if (nsites) {
      params.push(`n_sites=${nsites}`)
    }
    if (formula) {
      params.push(`chemical_formula_reduced=${formula}`)
    }
    resourcesApi.get(`/?${params.join('&')}`).then(response => {
        setExternalResources(response)
        setPagination(pagination => ({...pagination, total: response.data.length}))
      }).catch(raiseError)
    }, [archive, resourcesApi, raiseError, setExternalResources, setPagination])

  const columns = useMemo(() => [
    {
      key: 'id',
      sortable: true,
      align: 'left',
      render: data => (
        <Typography noWrap>
          <Ellipsis>{data.id}</Ellipsis>
        </Typography>
      )
    },
    {
      key: 'url',
      sortable: true,
      align: 'left',
      render: data => (
        <Typography noWrap>
          <Link href={data.url}>
            <Ellipsis>{data.url}</Ellipsis>
          </Link>
        </Typography>
      )
    },
    {
      key: 'database',
      sortable: true,
      align: 'left',
      render: data => data.database_name
    }
  ], [])

  if (!externalResources || (!externalResources.is_retrieving_more && externalResources.data.length === 0)) {
    return ''
  }

  const {page_size, page, total, order, order_by} = pagination
  const compare = (a, b) => {
    const value = data => String(data[order_by])
    return value(a).localeCompare(value(b)) * (order === 'asc' ? 1 : -1)
  }
  const start = page_size * (page - 1)
  const end = page_size + start
  const length = (end > total) ? total - start : page_size
  const data = [...externalResources.data].sort(compare).slice(start, start + length)

  return <PropertyCard title="Related resources">
    {externalResources.is_retrieving_more && externalResources.data.length === 0 && (
      <Box padding={2}>
        <Typography>Updating external resources, try again later</Typography>
      </Box>
    )}
    {(externalResources.data.length > 0) && (
      <Datatable
        columns={columns}
        data={data}
        pagination={pagination} onPaginationChanged={setPagination}
      >
        <DatatableTable details={ResourceDetails} actions={ResourceActions}>
          <DatatablePagePagination pageSizeValues={[10, 50, 100]}/>
        </DatatableTable>
      </Datatable>
    )}
  </PropertyCard>
})

RelatedResourcesCard.propTypes = {
  index: PropTypes.object.isRequired,
  archive: PropTypes.object
}

export default RelatedResourcesCard
