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
import React, {useMemo, useContext} from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import {
  makeStyles,
  Typography,
  Tooltip,
  IconButton,
  TableBody,
  TableContainer,
  TableHead,
  Table,
  TableRow,
  TableCell,
  Link,
  Box
} from '@material-ui/core'
import { DOI } from './dataset/DOI'
import ClipboardIcon from '@material-ui/icons/Assignment'
import { CopyToClipboard } from 'react-copy-to-clipboard'
import { get, isNil } from 'lodash'
import searchQuantities from '../searchQuantities'
import Placeholder from './visualization/Placeholder'
import NoData from './visualization/NoData'
import { formatNumber, formatTimestamp, authorList, serializeMetainfo } from '../utils'
import { Quantity as Q, Unit, useUnits } from '../units'
import { filterData } from './search/FilterRegistry'
import { RouteLink } from './nav/Routes'

/**
 * Component for showing a metainfo quantity value together with a name and
 * description.
*/
const useQuantityStyles = makeStyles(theme => ({
  root: {
    maxWidth: theme.spacing(35)
  },
  valueContainer: {
    display: 'flex',
    alignItems: 'center',
    flexDirection: 'row',
    maxWidth: '100%'
  },
  value: {
    flexGrow: 1
  },
  ellipsis: {
    whiteSpace: 'nowrap',
    overflow: 'hidden',
    textOverflow: 'ellipsis'
  },
  ellipsisFront: {
    direction: 'rtl',
    textAlign: 'left'
  },
  valueAction: {},
  valueActionButton: {
    padding: 4
  },
  valueActionIcon: {
    fontSize: 16
  },
  row: {
    display: 'flex',
    flexWrap: 'wrap',
    flexDirection: 'row',
    '& > :not(:last-child)': {
      marginRight: theme.spacing(3)
    }
  },
  column: {
    display: 'flex',
    flexDirection: 'column',
    '& > :not(:first-child)': {
      marginTop: theme.spacing(1)
    }
  },
  flex: {
    display: 'flex',
    flexDirection: 'row',
    flexWrap: 'wrap',
    alignContent: 'flex-start',
    '& div': {
      marginRight: theme.spacing(1)
    }
  },
  label: {
    color: 'rgba(0, 0, 0, 0.54)',
    fontSize: '0.75rem',
    fontWeight: 500
  },
  quantityList: {
    display: 'flex',
    flexDirection: 'column'
  }
}))

const Quantity = React.memo((props) => {
  const styles = useQuantityStyles()
  const presets = quantityPresets[props.quantity] || {}
  const {
    quantity,
    label,
    description,
    value,
    loading,
    placeholder,
    typography,
    noWrap,
    noLabel,
    row,
    column,
    flex,
    data,
    withClipboard,
    ellipsisFront,
    hideIfUnavailable,
    format
  } = {...presets, ...props}
  const children = props.children || (presets.render && presets.render(data))
  const units = useUnits()
  let content = null
  let clipboardContent = null

  let valueClassName = styles.value
  if (noWrap && ellipsisFront) {
    valueClassName = `${valueClassName} ${styles.ellipsisFront}`
  }

  const def = typeof quantity === 'string'
    ? searchQuantities[quantity]
    : undefined

  // Determine the final label to show
  const useLabel = useMemo(() => {
    let useLabel = label
    if (!useLabel) {
      if (def?.name) {
        useLabel = def.name.replace(/_/g, ' ')
      } else if (typeof quantity === 'string') {
        useLabel = quantity
      } else {
        useLabel = 'MISSING LABEL'
      }
    }
    return useLabel
  }, [quantity, label, def])

  const tooltip = description || def?.description || ''

  // Determine the final value to show.
  if (!loading) {
    let finalValue = value
    if (isNil(value)) {
      if (typeof quantity === 'string') {
        finalValue = data && quantity && get(data, quantity)
      } else if (children) {
      } else {
        try {
          finalValue = quantity(data)
        } catch {
          finalValue = undefined
        }
      }
    }

    if (finalValue === 'not processed') {
      finalValue = 'unavailable'
    }

    if (finalValue === 'unavailable') {
      finalValue = ''
    }

    if (format) {
      finalValue = serializeMetainfo(quantity, finalValue, units)
    }

    if ((!finalValue && !children) && hideIfUnavailable) {
      return null
    }

    if (finalValue && Array.isArray(finalValue)) {
      finalValue = finalValue.join(', ')
    }
    clipboardContent = clipboardContent || finalValue

    if (children && children.length !== 0) {
      content = children
    } else if (finalValue || finalValue === 0) {
      content = <Typography noWrap={noWrap} variant={typography} className={valueClassName}>
        {finalValue}
      </Typography>
    } else {
      content = <Typography noWrap={noWrap} variant={typography} className={valueClassName}>
        <i>{placeholder || 'unavailable'}</i>
      </Typography>
    }
  }

  if (row || column || flex) {
    return <div className={row ? styles.row : (column ? styles.column : styles.flex)}>{children}</div>
  } else {
    return (
      <div className={styles.root}>
        {!noLabel ? <Typography
          noWrap
          classes={{root: styles.label}}
          variant="caption"
        >{useLabel}</Typography> : ''}
        <div className={styles.valueContainer}>
          {loading
            ? <Typography noWrap={noWrap} variant={typography} className={valueClassName}>
              <i>loading ...</i>
            </Typography>
            // The tooltip portal is disabled for custom contents that may
            // contain links. Pressing a link while the tooltip is shown will
            // cause a navigation that leaves the popup opened (the Tooltip
            // state does not get updated since the page may be cached and a new
            // page is shown immediately).
            : <Tooltip title={tooltip} PopperProps={children ? {disablePortal: true} : undefined}>
              {content}
            </Tooltip>
          }
          {withClipboard
            ? <CopyToClipboard
              className={styles.valueAction}
              text={clipboardContent}
              onCopy={() => null}
            >
              <Tooltip title={`Copy ${useLabel} to clipboard`}>
                <div>
                  <IconButton
                    disabled={!clipboardContent}
                    classes={{root: styles.valueActionButton}}
                  >
                    <ClipboardIcon classes={{root: styles.valueActionIcon}}/>
                  </IconButton>
                </div>
              </Tooltip>
            </CopyToClipboard>
            : ''
          }
        </div>
      </div>
    )
  }
})

Quantity.propTypes = {
  children: PropTypes.node,
  label: PropTypes.string,
  typography: PropTypes.string,
  loading: PropTypes.bool,
  placeholder: PropTypes.string,
  noWrap: PropTypes.bool,
  noLabel: PropTypes.bool,
  row: PropTypes.bool,
  column: PropTypes.bool,
  flex: PropTypes.bool,
  data: PropTypes.object,
  value: PropTypes.any,
  quantity: PropTypes.oneOfType([
    PropTypes.string,
    PropTypes.func
  ]),
  withClipboard: PropTypes.bool,
  ellipsisFront: PropTypes.bool,
  hideIfUnavailable: PropTypes.bool,
  description: PropTypes.string,
  format: PropTypes.bool
}

export default Quantity

// Preset configuration for quantities
const quantityPresets = {
  datasets: {
    label: 'datasets',
    placeholder: 'no datasets',
    render: (data) => (data.datasets && data.datasets.length !== 0) &&
      <div>
        {data.datasets.map(ds => (
          <Typography key={ds.dataset_id}>
            <RouteLink path={`dataset/id/${ds.dataset_id}`}>{ds.dataset_name}</RouteLink>
            {ds.doi ? <span>&nbsp;<DOI style={{display: 'inline'}} parentheses doi={ds.doi}/></span> : ''}
          </Typography>))}
      </div>
  },
  authors: {
    label: 'authors',
    placeholder: 'no authors',
    render: (data) => <Typography>
      {authorList(data || [])}
    </Typography>
  },
  references: {
    label: 'references',
    placeholder: 'no references',
    render: (data) => (data?.references && <div style={{display: 'inline-grid'}}>
      {data.references.map(ref => <Typography key={ref} noWrap>
        <Link href={ref}>{ref}</Link>
      </Typography>)}
    </div>)
  },
  comment: {
    label: 'comment',
    placeholder: 'no comment'
  },
  upload_id: {
    noWrap: true,
    withClipboard: true,
    render: (data) => (
      <Box flexGrow={1}>
        <Typography noWrap>
          <RouteLink path={`upload/id/${data.upload_id}`}>{data.upload_id}</RouteLink>
        </Typography>
      </Box>
    )
  },
  last_processing_time: {
    noWrap: true,
    placeholder: 'not processed',
    render: (data) => <Typography noWrap>
      {formatTimestamp(data.last_processing_time)}
    </Typography>
  },
  last_processing_version: {
    description: 'Version used in the last processing',
    label: 'processing version',
    noWrap: true,
    placeholder: 'not processed',
    render: (data) => <Typography noWrap>
      {data.nomad_version}/{data.nomad_commit}
    </Typography>
  },
  upload_create_time: {
    noWrap: true,
    render: (data) => <Typography noWrap>
      {formatTimestamp(data.upload_create_time)}
    </Typography>
  },
  entry_id: {
    noWrap: true,
    withClipboard: true,
    render: (data) => (
      <Box flexGrow={1}>
        <Typography noWrap>
          <RouteLink path={`entry/id/${data.upload_id}/${data.entry_id}`}>{data.entry_id}</RouteLink>
        </Typography>
      </Box>
    )
  },
  'results.material.material_id': {
    noWrap: true,
    withClipboard: true
  },
  mainfile: {
    noWrap: true,
    ellipsisFront: true,
    withClipboard: true
  },
  raw_id: {
    noWrap: true,
    withClipboard: true,
    hideIfUnavailable: true
  },
  external_id: {
    noWrap: true,
    withClipboard: true,
    hideIfUnavailable: true
  }
}
/**
 * Representational component for tables containing metainfo data.
 */
const useTableStyles = makeStyles(theme => ({
  root: {
    border: `1px solid ${theme.palette.grey[300]}`
  }
}))
export const MetaInfoTable = React.memo(({data, className, classes, children}) => {
  const styles = useTableStyles(classes)
  return <TableContainer className={clsx(className, styles.root)}>
    <Table size="small">
      {children}
    </Table>
  </TableContainer>
})
MetaInfoTable.propTypes = {
  data: PropTypes.object,
  className: PropTypes.string,
  classes: PropTypes.object,
  children: PropTypes.node
}

/**
 * Used to organize individual quantities in a table.
 */
const quantityTableContext = React.createContext()
export const QuantityTable = React.memo(({data, className, children}) => {
  return <quantityTableContext.Provider value={data}>
    <MetaInfoTable className={className}>
      <TableBody>
        {children}
      </TableBody>
    </MetaInfoTable>
  </quantityTableContext.Provider>
})
QuantityTable.propTypes = {
  data: PropTypes.object,
  className: PropTypes.string,
  children: PropTypes.node
}

/**
 * Used to organize Quantities in a table row.
 */
const useRowStyles = makeStyles(theme => ({
  root: {}
}))
export const QuantityRow = React.memo(({className, classes, children}) => {
  const styles = useRowStyles()

  return <TableRow className={clsx(className, styles.root)}>
    {children}
  </TableRow>
})

QuantityRow.propTypes = {
  className: PropTypes.string,
  classes: PropTypes.object,
  children: PropTypes.node
}

/**
 * Used to display a quantity in a table cell.
 */
export const QuantityCell = React.memo(({
  quantity,
  value,
  data,
  label,
  description,
  classes,
  className,
  children,
  ...other
}) => {
  const contextData = useContext(quantityTableContext)
  const finalData = data || contextData

  return <TableCell align="left" {...other}>
    {children || <Quantity
      quantity={quantity}
      value={value}
      label={label}
      description={description}
      format
      noWrap
      data={finalData}
    />}
  </TableCell>
})

QuantityCell.propTypes = {
  quantity: PropTypes.oneOfType([PropTypes.string, PropTypes.func]),
  value: PropTypes.any,
  data: PropTypes.object,
  label: PropTypes.string,
  description: PropTypes.string,
  options: PropTypes.object,
  className: PropTypes.string,
  classes: PropTypes.object,
  children: PropTypes.node
}

/**
 * Used to display data from one or many sections in a table.
 */
const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    height: '100%'
  },
  table: {
    marginBottom: theme.spacing(1)
  }
}))
export const SectionTable = React.memo(({
  data,
  section,
  quantities,
  horizontal,
  showIndex,
  classes,
  className,
  units,
  'data-testid': testID
}) => {
  const styles = useStyles({classes: classes})
  if (showIndex && !horizontal) {
    throw Error('Showing index in a vertically displayed table is not currently supported.')
  }

  // If data is set explicitly to False, we show the NoData component.
  let content
  if (data === false) {
    content = <NoData data-testid={`${testID}-nodata`}/>
  } else if (!data) {
    content = <Placeholder variant="rect" data-testid={`${testID}-placeholder`}/>
  } else {
    content = <MetaInfoTable className={styles.table}>
      <TableHead>
        {horizontal
          ? <TableRow>
            {showIndex && <TableCell align="left">
              <Tooltip title="Item index">
                <span>Index</span>
              </Tooltip>
            </TableCell>}
            {Object.keys(quantities).map((key, index) => {
              const defCustom = quantities[key]
              const def = filterData[`${section}.${key}`]
              const unitName = defCustom.unit || def?.unit
              const unit = unitName && new Unit(unitName)
              const unitLabel = unit && unit.toSystem(units).label()
              const label = defCustom.label || def?.label
              const description = defCustom.description || def?.description || ''
              const content = unit ? `${label} (${unitLabel})` : defCustom.label
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
              {showIndex && <TableCell align="left">{i}</TableCell>}
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
                      unit ? new Q(value, unit).toSystem(units).value() : value,
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
          : <>{Object.keys(quantities).map((key, i) => (
            <TableRow key={i}>
              {data.data.map((row, j) => {
                const defCustom = quantities[key]
                const def = filterData[`${section}.${key}`]
                const unitName = defCustom.unit || def?.unit
                const unit = unitName && new Unit(unitName)
                const unitLabel = unit ? ` ${unit.toSystem(units).label()}` : ''
                const description = defCustom.description || def?.description || ''
                const dtype = defCustom?.type?.type_data || def?.type?.type_data
                const align = defCustom.align || 'right'
                const label = defCustom.label || def?.label
                let value = row[key]
                if (value !== undefined) {
                  if (!isNaN(value)) {
                    value = `${formatNumber(
                      unit ? new Q(value, unit).toSystem(units).value() : value,
                      dtype
                    )}${unitLabel}`
                  }
                } else {
                  value = defCustom.placeholder || 'unavailable'
                }
                return <>
                  <TableCell key={j} align={align}>
                    <Tooltip title={description}>
                      <span>
                        {label}
                      </span>
                    </Tooltip>
                  </TableCell>
                  <TableCell key={j} align={align}>{value}</TableCell>
                </>
              })}
            </TableRow>
          ))}</>
        }
      </TableBody>
    </MetaInfoTable>
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
  showIndex: PropTypes.bool, // Whether to display the item index
  className: PropTypes.string,
  classes: PropTypes.object,
  units: PropTypes.object,
  'data-testid': PropTypes.string
}
