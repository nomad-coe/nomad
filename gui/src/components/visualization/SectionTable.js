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
import {
  Tooltip,
  Table,
  TableBody,
  TableRow,
  TableHead,
  TableContainer,
  TableCell
} from '@material-ui/core'
import { makeStyles } from '@material-ui/core/styles'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import Placeholder from './Placeholder'
import NoData from './NoData'
import { formatNumber } from '../../utils'
import { Unit, toUnitSystem } from '../../units'
import searchQuantities from '../../searchQuantities'

/**
 * Used to display data from one or many sections in a table.
 */
const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    height: '100%'
  },
  table: {
    width: '100%',
    height: '100%',
    marginBottom: theme.spacing(1)
  }
}))
const SectionTable = React.memo(({
  data,
  section,
  quantities,
  horizontal,
  classes,
  className,
  units,
  'data-testid': testID
}) => {
  const styles = useStyles({classes: classes})

  // If data is set explicitly to False, we show the NoData component.
  let content
  if (data === false) {
    content = <NoData data-testid={`${testID}-nodata`}/>
  } else if (!data) {
    content = <Placeholder variant="rect" data-testid={`${testID}-placeholder`}/>
  } else {
    content = <TableContainer className={styles.table}>
      <Table size="small" aria-label="table">
        <TableHead>
          {horizontal
            ? <TableRow>
              {Object.keys(quantities).map((key, index) => {
                const defCustom = quantities[key]
                const def = searchQuantities[`${section}.${key}`]
                const unitName = defCustom.unit || def?.unit
                const unit = unitName && new Unit(unitName)
                const unitLabel = unit && unit.label(units)
                const description = defCustom.description || def.description || ''
                const content = unit ? `${defCustom.label} (${unitLabel})` : defCustom.label
                const align = defCustom.align || 'right'
                return <TableCell key={index} align={align}>
                  <Tooltip title={description}>
                    <span>
                      {content}
                    </span>
                  </Tooltip>
                </TableCell>
              })}
            </TableRow>
            : null
          }
        </TableHead>
        <TableBody>
          {horizontal
            ? <>{data.data.map((row, i) => (
              <TableRow key={i}>
                {Object.keys(quantities).map((key, j) => {
                  const defCustom = quantities[key]
                  const def = searchQuantities[`${section}.${key}`]
                  const unit = defCustom.unit || def?.unit
                  const dtype = defCustom?.type?.type_data || def?.type?.type_data
                  const align = defCustom.align || 'right'
                  let value = row[key]
                  if (value !== undefined) {
                    if (!isNaN(value)) {
                      value = formatNumber(
                        unit ? toUnitSystem(value, unit, units, false) : value,
                        dtype
                      )
                    }
                  } else {
                    value = defCustom.placeholder || 'unavailable'
                  }
                  return <TableCell key={j} align={align}>{value}</TableCell>
                })}
              </TableRow>
            ))}</>
            : null
          }
        </TableBody>
      </Table>
    </TableContainer>
  }

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    {content}
  </div>
})

SectionTable.propTypes = {
  data: PropTypes.oneOfType([
    PropTypes.bool, // Set to False to show NoData component
    PropTypes.shape({
      data: PropTypes.arrayOf(PropTypes.object).isRequired
    })
  ]),
  section: PropTypes.string,
  quantities: PropTypes.any,
  horizontal: PropTypes.bool,
  className: PropTypes.string,
  classes: PropTypes.object,
  units: PropTypes.object,
  'data-testid': PropTypes.string
}

export default SectionTable
