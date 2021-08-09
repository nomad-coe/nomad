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
  Typography,
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
import AspectRatio from './AspectRatio'
import { formatNumber } from '../../utils'
import { Unit, toUnitSystem } from '../../units'
import searchQuantities from '../../searchQuantities'

/**
 * Used to display data from one or many sections in a table.
 */
const useStyles = makeStyles(theme => ({
  root: {},
  table: {
    padding: '1em',
    boxSizing: 'border-box'
  }
}))
const SectionTable = React.memo(({
  data,
  label,
  section,
  quantities,
  horizontal,
  aspectRatio,
  classes,
  className,
  units,
  'data-testid': testID
}) => {
  const styles = useStyles({classes: classes})
  return <AspectRatio
    aspectRatio={aspectRatio}
    className={clsx(className, styles.root)}
    data-testid={testID}
  >
    { data
      ? <TableContainer className={styles.table}>
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
                  return <TableCell key={index} align="left">
                    <Tooltip title={description}>
                      <Typography variant="overline">{content}</Typography>
                    </Tooltip>
                  </TableCell>
                })}
              </TableRow>
              : null
            }
          </TableHead>
          <TableBody>
            {horizontal
              ? <>{data.map((row, i) => (
                <TableRow key={i}>
                  {Object.keys(quantities).map((key, j) => {
                    const defCustom = quantities[key]
                    const def = searchQuantities[`${section}.${key}`]
                    const unit = defCustom.unit || def?.unit
                    const dtype = defCustom?.type?.type_data || def?.type?.type_data
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
                    return <TableCell key={j} align="left">{value}</TableCell>
                  })}
                </TableRow>
              ))}</>
              : null
            }
          </TableBody>
        </Table>
      </TableContainer>
      : <Placeholder variant="rect" data-testid={`${testID}-placeholder`}/>
    }
  </AspectRatio>
})

SectionTable.propTypes = {
  data: PropTypes.any,
  label: PropTypes.string,
  section: PropTypes.string,
  quantities: PropTypes.any,
  horizontal: PropTypes.bool,
  aspectRatio: PropTypes.number,
  className: PropTypes.string,
  classes: PropTypes.object,
  units: PropTypes.object,
  'data-testid': PropTypes.string
}

export default SectionTable
