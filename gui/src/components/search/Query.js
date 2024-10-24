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
import React, { useMemo } from 'react'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { isNil, isEmpty, isPlainObject } from 'lodash'
import { useSearchContext } from './SearchContext'
import { useUnitContext } from '../units/UnitContext'
import { Typography, Box, Chip, Tooltip } from '@material-ui/core'
import FilterTitle from './FilterTitle'
import Ellipsis from '../visualization/Ellipsis'
import ClearIcon from '@material-ui/icons/Clear'
import ReplayIcon from '@material-ui/icons/Replay'
import LockIcon from '@material-ui/icons/Lock'
import CodeIcon from '@material-ui/icons/Code'
import { Actions, Action } from '../Actions'
import { SourceApiCall, SourceApiDialogButton, SourceDialogDivider, SourceJsonCode } from '../buttons/SourceDialogButton'

/**
 * Thin wrapper for MUI Chip that is used for displaying (and possibly removing)
 * query values.
 */
const useQueryChipStyles = makeStyles(theme => ({
  root: {
    maxWidth: '100%'
  },
  chipRoot: {
    width: '100%',
    maxWidth: '100%'
  },
  chipLabel: {
    minWidth: '1rem',
    maxWidth: '25rem'
  }
}))
export const QueryChip = React.memo(({
  label,
  onDelete,
  color,
  className,
  locked
}) => {
  const styles = useQueryChipStyles()

  return <div className={clsx(className, styles.root)}>
    <Chip
      label={<Ellipsis tooltip={label}>{label}</Ellipsis>}
      onDelete={locked ? undefined : onDelete}
      color={locked ? undefined : color}
      icon={locked ? <LockIcon/> : undefined}
      classes={{root: styles.chipRoot, label: styles.chipLabel}}
    />
  </div>
})

QueryChip.propTypes = {
  label: PropTypes.oneOfType([PropTypes.string, PropTypes.number]),
  onDelete: PropTypes.func,
  color: PropTypes.string,
  className: PropTypes.string,
  locked: PropTypes.bool
}
QueryChip.defaultProps = {
  color: 'primary'
}

/**
 * Used to group several related query chips inside one container.
 */
export const queryTitleHeight = 2.2
export const queryGroupHeight = 4.1 + queryTitleHeight
export const queryGroupSpacing = 0.5
const useQueryChipGroupStyles = makeStyles(theme => ({
  root: {
    position: 'relative',
    minHeight: theme.spacing(queryGroupHeight),
    marginLeft: theme.spacing(queryGroupSpacing),
    marginRight: theme.spacing(queryGroupSpacing)
  },
  chips: {
    marginTop: theme.spacing(queryTitleHeight),
    display: 'flex',
    flexDirection: 'row',
    flexWrap: 'wrap',
    alignItems: 'center',
    justifyContent: 'center',
    backgroundColor: theme.palette.primary.main,
    borderRadius: theme.spacing(2)
  },
  titleRoot: {
    position: 'absolute',
    left: theme.spacing(0.4),
    right: 0,
    top: 0,
    height: theme.spacing(queryTitleHeight)
  },
  title: {
    color: theme.palette.grey[600]
  }
}))
export const QueryChipGroup = React.memo(({
  quantity,
  className,
  children
}) => {
  const styles = useQueryChipGroupStyles()
  return <div className={clsx(className, styles.root)}>
    <FilterTitle
      quantity={quantity}
      variant="caption"
      classes={{root: styles.titleRoot, text: styles.title}}
      disableUnit
    />
    <div className={styles.chips}>
      {children}
    </div>
  </div>
})

QueryChipGroup.propTypes = {
  quantity: PropTypes.string,
  color: PropTypes.string,
  className: PropTypes.string,
  children: PropTypes.node
}
QueryChipGroup.defaultProps = {
  color: 'primary'
}

/**
 * Operators between query chips.
 */
const useQueryOpStyles = makeStyles(theme => ({
  root: {
    cursor: 'default'
  },
  inside: {
    fontSize: '0.6rem',
    marginLeft: theme.spacing(0.2),
    marginRight: theme.spacing(0),
    color: 'white'
  },
  outside: {
    fontSize: '1.5rem',
    marginLeft: theme.spacing(0.1),
    marginRight: theme.spacing(0.1),
    marginBottom: theme.spacing(-0.4)
  }
}))
export const QueryOp = React.memo(({className, children, tooltip, variant}) => {
  const styles = useQueryOpStyles()
  return <Tooltip title={tooltip || ''}>
    <Typography
      variant="caption"
      className={clsx(className, styles.root, styles[variant])}
    >{children}
    </Typography>
  </Tooltip>
})

QueryOp.propTypes = {
  className: PropTypes.string,
  tooltip: PropTypes.string,
  children: PropTypes.node,
  variant: PropTypes.string
}
QueryOp.defaultProps = {
  variant: 'inside'
}

export const QueryAnd = React.memo(() => {
  return <QueryOp>AND</QueryOp>
})
export const QueryOr = React.memo(() => {
  return <QueryOp>OR</QueryOp>
})
export const QueryCurlyBracketLeft = React.memo((props) => {
  return <QueryOp tooltip="Starts a nested query" {...props}>{'{'}</QueryOp>
})
export const QueryCurlyBracketRight = React.memo((props) => {
  return <QueryOp tooltip="Ends a nested query" {...props}>{'}'}</QueryOp>
})

// Custom function for chip creation
const createChips = (name, filterValue, onDelete, filterData, units) => {
  if (isNil(filterValue)) return []

  const { serializerPretty: serializer, customSerialization, queryMode } = filterData[name]
  const isArray = Array.isArray(filterValue)
  const isSet = filterValue instanceof Set
  const isObj = isPlainObject(filterValue)
  const op = queryMode === "any" ? <QueryOr/> : <QueryAnd/>
  const chips = []

  const createChip = (label, onDelete, single = false) => (
    <QueryChip key={label} label={label} onDelete={onDelete} single={single}/>
  )

  if (customSerialization) {
    chips.push({ comp: createChip(serializer(filterValue), () => onDelete(undefined)), op })
  } else if (isArray || isSet) {
    Array.from(filterValue).forEach((value, index) => {
      chips.push({
        comp: createChip(serializer(value, units), () => {
          const newValue = isSet ? new Set(filterValue) : [...filterValue]
          isSet ? newValue.delete(value) : newValue.splice(index, 1)
          onDelete(newValue)
        }),
        op
      })
    })
  } else if (isObj) {
    const createRangeChip = (label, comparison) => {
      const content = `${comparison} ${serializer(filterValue[label], units)}`
      if (!isNil(filterValue[label])) {
        let newValue = {...filterValue}
        delete newValue[label]
        newValue = isEmpty(newValue) ? undefined : newValue
        chips.push({ comp: createChip(content, () => onDelete(newValue)), op })
      }
    }

    createRangeChip('gte', '>=')
    createRangeChip('gt', '>')
    createRangeChip('lte', '<=')
    createRangeChip('lt', '<')
  } else {
    chips.push({ comp: createChip(serializer(filterValue), () => onDelete(undefined)), op })
  }

  return chips.length ? (
    <QueryChipGroup key={name} quantity={name}>
      {chips.map((chip, index) => (
        <React.Fragment key={index}>
          {chip.comp}
          {index < chips.length - 1 && chip.op}
        </React.Fragment>
      ))}
    </QueryChipGroup>
  ) : null
}

/*
 * Displays chips for the current query.
 */
const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    boxSizing: 'border-box',
    display: 'flex',
    flexDirection: 'row',
    flexWrap: 'wrap',
    alignItems: 'flex-end',
    marginLeft: theme.spacing(-queryGroupSpacing),
    marginRight: theme.spacing(-queryGroupSpacing)
  },
  empty: {
    marginTop: theme.spacing(1.8)
  },
  chip: {
    padding: theme.spacing(0.5)
  }
}))
const QueryChips = React.memo(({ className, classes }) => {
  const { filterData, useQuery, useUpdateFilter } = useSearchContext()
  const query = useQuery()
  const updateFilter = useUpdateFilter()
  const theme = useTheme()
  const { units } = useUnitContext()
  const styles = useStyles({ classes, theme })

  const chips = useMemo(() => {
    const chips = []
    // The query chips are created in alphabetical order
    const keys = Object.keys(query)
    keys.sort()
    for (const quantity of keys) {
      const filterValue = query[quantity]
      // Each key in a section is mapped into a group
      const isSection = filterData[quantity].section
      if (isSection) {
        const addChipsForSection = (data) => {
          const newChips = []
          Object.entries(data).forEach(([key, value], index) => {
            // Empty filters are skipped
            if (isEmpty(value)) return
            const onDelete = (newValue) => {
              const newSection = { ...data, [key]: newValue }
              if (newValue === undefined) delete newSection[key]
              updateFilter([quantity, newSection])
            }
            newChips.push(
              <React.Fragment key={`${quantity}.${key}`}>
                {createChips(`${quantity}.${key}`, value, onDelete, filterData, units)}
              </React.Fragment>
            )
          })
          return newChips
        }
        const newChips = addChipsForSection(filterValue)
        if (newChips.length) {
          chips.push(<QueryCurlyBracketLeft variant='outside'/>)
          chips.push(newChips)
          chips.push(<QueryCurlyBracketRight variant='outside'/>)
        }
      // Regular chips get their own group
      } else {
        const onDelete = (newValue) => updateFilter([quantity, newValue])
        chips.push(createChips(quantity, filterValue, onDelete, filterData, units))
      }
    }

    return chips
  }, [query, filterData, units, updateFilter])

  return (
    <div className={clsx(className, styles.root)}>
      {chips.length
        ? chips
        : <Typography className={styles.empty}><i>Your query will be shown here</i></Typography>
      }
    </div>
  )
})

QueryChips.propTypes = {
  className: PropTypes.string,
  classes: PropTypes.object
}

const useQueryStyles = makeStyles(theme => ({
  root: {
    minHeight: theme.spacing(queryGroupHeight),
    margin: theme.spacing(1, 0.25),
    width: '100%',
    display: 'flex',
    flexDirection: 'row',
    alignItems: 'flex-start'
  },
  offset: {
    height: theme.spacing(queryGroupHeight),
    display: 'flex',
    flexDirection: 'row',
    alignItems: 'center',
    width: 'unset'
  }
}))

/**
 * Displays the current query and actions for it.
 */
const Query = React.memo(() => {
  const {useResetFilters, useRefresh, useApiData} = useSearchContext()
  const styles = useQueryStyles()
  const resetFilters = useResetFilters()
  const refresh = useRefresh()
  const apiData = useApiData()

  return <Box className={styles.root}>
    <QueryChips/>
    <Actions className={styles.offset}>
      <Action
        tooltip="Clear query"
        onClick={() => resetFilters()}
      >
        <ClearIcon fontSize="small"/>
      </Action>
      <Action
        tooltip="Refresh results"
        onClick={() => refresh()}
      >
        <ReplayIcon fontSize="small"/>
      </Action>
      <Action
        tooltip=""
        ButtonComponent={SourceApiDialogButton}
        ButtonProps={{
          tooltip: "View API call for the query",
          maxWidth: "lg",
          fullWidth: true,
          icon: <CodeIcon fontSize="small"/>,
          buttonProps: {
            size: "small"
          }
        }}
      >
        <Typography>
          NOMAD uses the same query format throughout its API. This is the query
          based on the current filters:
        </Typography>
        <SourceJsonCode data={{owner: apiData?.body?.owner, query: apiData?.body?.query}}/>
        <SourceDialogDivider/>
        <Typography>
          One application of the above query is this API call. This is what is currently
          used to render this page and includes all displayed statistics data
          (aggregations).
        </Typography>
        <SourceApiCall
          {...apiData}
        />
      </Action>
    </Actions>
  </Box>
})
Query.propTypes = {}

export default Query
