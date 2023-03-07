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
import React, {useMemo, useContext, useCallback} from 'react'
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
import { get, isNil, isString, isNumber } from 'lodash'
import searchQuantities from '../searchQuantities'
import Placeholder from './visualization/Placeholder'
import NoData from './visualization/NoData'
import { formatNumber, formatTimestamp, authorList, serializeMetainfo } from '../utils'
import { Quantity as Q, Unit, useUnits } from '../units'
import { filterData } from './search/FilterRegistry'
import { MaterialLink, RouteLink } from './nav/Routes'

/**
 * Component for showing a metainfo quantity value together with a name and
 * description.
*/
const useQuantityStyles = makeStyles(theme => ({
  root: {
    width: '100%'
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
    maxWidth,
    format
  } = {...presets, ...props}

  const getRenderFromType = useCallback((quantity, data) => {
    const type = quantity.type
    if (type.type_data === 'str') {
      return <Typography noWrap>
        {data[quantity.name]}
      </Typography>
    } else if (type.type_data === 'nomad.metainfo.metainfo._Datetime') {
      return <Typography noWrap>
        {formatTimestamp(data[quantity.name])}
      </Typography>
    }
    return <Typography noWrap>
      {data[quantity.name]}
    </Typography>
  }, [])

  // Determine the rendered children. They may have been given explicitly, or
  // this quantity may have a default rendering function for a value, or this
  // quantity may have a default rendering function for data, or the data type
  // may be associated with a particular rendering function.
  const children = props.children ||
    ((presets.render && !isNil(data)) && presets.render(data)) ||
    ((presets.renderValue && !isNil(value)) && presets.renderValue(value)) ||
    ((quantity?.name && !isNil(quantity.type)) && getRenderFromType(quantity, data))

  const units = useUnits()
  let content = null
  let clipboardContent = null

  let valueClassName = styles.value
  if (noWrap && ellipsisFront) {
    valueClassName = `${valueClassName} ${styles.ellipsisFront}`
  }

  const isQuantityString = isString(quantity)
  const def = isQuantityString
    ? (searchQuantities[quantity])
    : quantity?.name && quantity.type && quantity

  // Determine the final label to show
  const useLabel = useMemo(() => {
    let useLabel = label
    if (!useLabel) {
      // Primarily use a lowercase 'pretty' label if one is defined in FilterRegistry
      if (isQuantityString && filterData?.[quantity]?.label) {
        useLabel = filterData?.[quantity]?.label.toLowerCase()
      // Alternatively use the original name in metainfo, underscores replaced by spaces
      } else if (def?.name) {
        useLabel = def.name.replace(/_/g, ' ')
      } else if (isQuantityString) {
        useLabel = quantity
      } else {
        useLabel = 'MISSING LABEL'
      }
    }
    return useLabel
  }, [quantity, label, def, isQuantityString])

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
    return <div className={row ? styles.row : (column ? styles.column : styles.flex)} data-testid={`quantity-${def?.name}`}>{children}</div>
  } else {
    return (
      <div className={styles.root} style={{maxWidth: maxWidth}} data-testid={noLabel ? undefined : `quantity-${def?.name}`}>
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
    PropTypes.func,
    PropTypes.object
  ]),
  withClipboard: PropTypes.bool,
  ellipsisFront: PropTypes.bool,
  hideIfUnavailable: PropTypes.bool,
  description: PropTypes.string,
  format: PropTypes.bool,
  maxWidth: PropTypes.string
}

Quantity.defaultProps = {
  maxWidth: '280px'
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
          <RouteLink path={`entry/id/${data.entry_id}`}>{data.entry_id}</RouteLink>
        </Typography>
      </Box>
    )
  },
  'results.material.material_id': {
    noWrap: true,
    hideIfUnavailable: true,
    placeholder: 'unavailable',
    withClipboard: true,
    render: (data) => {
      return data?.results?.material?.material_id
        ? <Box flexGrow={1}>
          <Typography noWrap>
            <MaterialLink materialId={data.results.material.material_id}>{data.results.material.material_id}</MaterialLink>
          </Typography>
        </Box>
        : null
    }
  },
  'results.material.topology.material_id': {
    noWrap: true,
    hideIfUnavailable: true,
    placeholder: 'unavailable',
    withClipboard: true,
    renderValue: (value) => {
      return value
        ? <Box flexGrow={1}>
          <Typography noWrap>
            <MaterialLink materialId={value}>{value}</MaterialLink>
          </Typography>
        </Box>
        : null
    }
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
const useMetaInfoTableStyles = makeStyles(theme => ({
  root: {
    border: `1px solid ${theme.palette.grey[300]}`
  },
  wrap: {
    borderBottom: `none`
  },
  fixed: {
    tableLayout: 'fixed'
  }
}))
export const MetaInfoTable = React.memo(({fixed, wrap, className, classes, children}) => {
  const styles = useMetaInfoTableStyles(classes)
  return <TableContainer className={clsx(className, styles.root, wrap && styles.wrap)}>
    <Table size="small" className={clsx(fixed && styles.fixed)}>
      {children}
    </Table>
  </TableContainer>
})
MetaInfoTable.propTypes = {
  fixed: PropTypes.bool, // Whether the table cells should have equal sizes
  wrap: PropTypes.bool, // Whether the table rows should wrap
  className: PropTypes.string,
  classes: PropTypes.object,
  children: PropTypes.node
}

/**
 * Used to organize individual quantities in a table.
 */
const quantityTableContext = React.createContext()
export const QuantityTable = React.memo(({data, fixed, wrap, className, children}) => {
  return <quantityTableContext.Provider value={{data, wrap}}>
    <MetaInfoTable wrap={wrap} fixed={fixed} className={className}>
      <TableBody>
        {children}
      </TableBody>
    </MetaInfoTable>
  </quantityTableContext.Provider>
})
QuantityTable.propTypes = {
  data: PropTypes.object,
  fixed: PropTypes.bool,
  wrap: PropTypes.bool,
  className: PropTypes.string,
  children: PropTypes.node
}

/**
 * Used to organize Quantities in a table row.
 */
const useRowStyles = makeStyles(theme => ({
  root: {},
  wrap: {
    display: 'inline-flex',
    flexWrap: 'wrap',
    overflow: 'hidden'
  }
}))
export const QuantityRow = React.memo(({className, classes, children}) => {
  const styles = useRowStyles()
  const wrap = useContext(quantityTableContext)?.wrap
  return <TableRow className={clsx(className, styles.root, wrap && styles.wrap)}>
    {children}
  </TableRow>
})

QuantityRow.propTypes = {
  className: PropTypes.string,
  wrap: PropTypes.bool, // Whether to wrap the cells once they overflow
  classes: PropTypes.object,
  children: PropTypes.node
}

/**
 * Used to display a quantity in a table cell.
 */
const useCellStyles = makeStyles(theme => ({
  wrap: {
    borderBottom: `1px solid ${theme.palette.grey[300]} !important`,
    flexGrow: 1
  }
}))
export const QuantityCell = React.memo(({
  quantity,
  value,
  data,
  label,
  description,
  classes,
  className,
  hideIfUnavailable,
  format,
  children,
  maxWidth,
  ...other
}) => {
  const context = useContext(quantityTableContext)
  const wrap = context?.wrap
  const contextData = context?.data
  const finalData = data || contextData
  const styles = useCellStyles()

  return (hideIfUnavailable && isNil(value))
    ? null
    : <TableCell align="left" {...other} className={clsx(wrap && styles.wrap)}>
      {children || <Quantity
        quantity={quantity}
        value={value}
        label={label}
        description={description}
        hideIfUnavailable={hideIfUnavailable}
        format={format}
        noWrap
        data={finalData}
        maxWidth={maxWidth}
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
  maxWidth: PropTypes.string,
  className: PropTypes.string,
  hideIfUnavailable: PropTypes.bool,
  format: PropTypes.bool,
  classes: PropTypes.object,
  children: PropTypes.node
}

QuantityCell.defaultProps = {
  format: true
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

/**
 * Recursively displays all of the quantities that are present within the given
 * section.
 */
export const SectionTableAutomatic = React.memo(({data, prefix, columns}) => {
  // Figure out a flat list of quantities to display
  const quantities = []
  const addMethodQuantities = (obj, parentKey) => {
    const children = {}
    Object.keys(obj).forEach(key => {
      const value = obj[key]
      if (Array.isArray(value) || isString(value) || isNumber(value)) {
        const metainfoPath = `${prefix}.${parentKey === '' ? '' : `${parentKey}.`}${key}`
        quantities.push({
          label: key.replace(/_/g, ' '),
          quantity: metainfoPath,
          value: value
        })
      } else if (value instanceof Object) {
        children[key] = value
      }
    })
    Object.keys(children).forEach(key => addMethodQuantities(children[key], `${parentKey === '' ? '' : `${parentKey}.`}${key}`))
  }
  addMethodQuantities(data, '')

  // Create rows that each have a fixed number of items
  const rows = []
  quantities.forEach((quantity, index) => {
    if (index % columns === 0) {
      rows.push([])
    }
    const row = rows[Math.floor(index / columns)]
    row.push(quantity)
  })

  return <QuantityTable>
    {rows.map((row, i) => <QuantityRow key={i}>
      {row.map((cell) => <QuantityCell
        key={cell.label}
        {...cell}
      />)}
    </QuantityRow>)}
  </QuantityTable>
})

SectionTableAutomatic.propTypes = {
  data: PropTypes.object,
  prefix: PropTypes.string, // Possible prefix that is needed to resolve the quantity metainfo
  columns: PropTypes.number // The number of columns to display
}

SectionTableAutomatic.defaultProps = {
  columns: 3
}
