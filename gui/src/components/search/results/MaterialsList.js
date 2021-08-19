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
import { makeStyles, TableCell, Toolbar, IconButton, Tooltip } from '@material-ui/core'
import NextIcon from '@material-ui/icons/ChevronRight'
import StartIcon from '@material-ui/icons/SkipPrevious'
import DataTable from '../../DataTable'
import DetailsIcon from '@material-ui/icons/MoreHoriz'
import { appBase } from '../../../config'

const useStyles = makeStyles(theme => ({
  root: {
    overflow: 'auto',
    paddingLeft: theme.spacing(2),
    paddingRight: theme.spacing(2)
  },
  scrollCell: {
    padding: 0
  },
  scrollBar: {
    minHeight: 56,
    padding: 0
  },
  scrollSpacer: {
    flexGrow: 1
  },
  clickableRow: {
    cursor: 'pointer'
  }
}))

const columns = {
  'encyclopedia.material.formula': {
    label: 'Formula'
  },
  'encyclopedia.material.material_name': {
    label: 'Name'
  },
  'encyclopedia.material.material_type': {
    label: 'Type'
  },
  'encyclopedia.material.bulk': {
    label: 'Spacegroup',
    render: entry => {
      const bulk = entry.encyclopedia.material.bulk
      return (bulk && bulk.space_group_international_short_symbol) || '-'
    }
  },
  calculations: {
    label: 'No calculations',
    description: 'The number of entries with data for this material',
    render: entry => entry.total
  }
}

export default function MaterialsList(props) {
  const { data, total, materials_after, per_page, onChange, actions } = props
  const classes = useStyles()
  const materials = data['encyclopedia.material.materials_grouped'] || {values: []}
  const results = Object.keys(materials.values).map(id => {
    return {
      id: id,
      total: materials.values[id].total,
      ...materials.values[id].examples[0]
    }
  })
  const after = materials.after
  const perPage = per_page || 10

  let paginationText
  if (materials_after) {
    paginationText = `next ${results.length.toLocaleString()} of ${(total || 0).toLocaleString()}`
  } else {
    paginationText = `1-${results.length.toLocaleString()} of ${(total || 0).toLocaleString()}`
  }

  /*   const handleViewMaterial = useCallback((event, materialId) => {
    event.stopPropagation()
    history.push(`/material/${materialId}/overview`)
  }, [history]) */

  const pagination = <TableCell colSpan={1000} classes={{root: classes.scrollCell}}>
    <Toolbar className={classes.scrollBar}>
      <span className={classes.scrollSpacer}>&nbsp;</span>
      <span>{paginationText}</span>
      <IconButton disabled={!materials_after} onClick={() => onChange({materials_grouped_after: null})}>
        <StartIcon />
      </IconButton>
      <IconButton disabled={results.length < perPage} onClick={() => onChange({materials_grouped_after: after})}>
        <NextIcon />
      </IconButton>
    </Toolbar>
  </TableCell>

  const entryActions = entry => <Tooltip title="Open this material in the Encyclopedia.">
    <IconButton href={`${appBase}/encyclopedia/#/material/${entry.encyclopedia.material.material_id}`}>
      <DetailsIcon />
    </IconButton>
    {/*     <IconButton onClick={event => handleViewMaterial(event, entry.encyclopedia.material.material_id)}>
      <DetailsIcon />
    </IconButton> */}
  </Tooltip>

  return <DataTable
    entityLabels={['material', 'materials']}
    id={row => row.id}
    total={total}
    columns={columns}
    selectedColumns={['encyclopedia.material.formula', 'encyclopedia.material.material_name', 'encyclopedia.material.material_type', 'encyclopedia.material.bulk.spacegroup', 'calculations']}
    selectedColumnsKey="materials"
    data={results}
    rows={perPage}
    actions={actions}
    pagination={pagination}
    entryActions={entryActions}
  />
}
MaterialsList.propTypes = ({
  data: PropTypes.object,
  total: PropTypes.number,
  onChange: PropTypes.func.isRequired,
  materials_after: PropTypes.string,
  per_page: PropTypes.number,
  actions: PropTypes.element
})
