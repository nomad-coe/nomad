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
import React, { useCallback, useState } from 'react'
import PropTypes from 'prop-types'
import { FilterSubMenu } from './FilterMenu'
import { InputGrid, InputGridItem } from '../input/InputGrid'
import { Checkbox, Collapse, makeStyles, Typography } from '@material-ui/core'
import { rootSections, resolveRef } from '../../archive/metainfo'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ChevronRightIcon from '@material-ui/icons/ChevronRight'
import { useSearchContext } from '../SearchContext'

const filterProperties = def => !(def.name.startsWith('x_') || def.virtual)

const useDefinitionStyles = makeStyles(theme => ({
  root: {},
  item: {
    '&:hover': {
      backgroundColor: theme.palette.grey[100]
    },
    width: '100%',
    display: 'flex',
    flexDirection: 'row',
    alignItems: 'center',
    flexWrap: 'nowrap',
    height: 32,
    paddingRight: theme.spacing(1)
  },
  icon: {
  },
  name: {
    flexGrow: 1
  },
  actions: {
    marginLeft: theme.spacing(1)
  },
  children: {
    marginLeft: theme.spacing(3)
  }
}))

const Definition = React.memo(function Definition({def, name, path}) {
  const classes = useDefinitionStyles()
  const [open, setOpen] = useState(false)
  const { useFilterState } = useSearchContext()
  const [filter, setFilter] = useFilterState('quantities')

  name = name || def.name

  const handleSelect = useCallback((event) => {
    event.stopPropagation()
    setFilter(old => {
      let newValues
      if (old) {
        const isSelected = old?.has(path)
        isSelected ? old.delete(path) : old.add(path)
        newValues = new Set(old)
      } else {
        newValues = new Set([path])
      }
      return newValues
    })
  }, [path, setFilter])

  const handleToggle = useCallback(event => {
    event.stopPropagation()
    setOpen(open => !open)
  }, [setOpen])

  const hasChildren = def.m_def === 'Section' && def.sub_sections.length + def.quantities.length

  const iconProps = {
    className: classes.icon,
    fontSize: 'small'
  }

  const icon = hasChildren
    ? (open ? <ExpandMoreIcon {...iconProps} /> : <ChevronRightIcon {...iconProps}/>)
    : <div className={classes.icon} />

  const childPathPrefix = path ? path + '.' : ''

  return <div className={classes.root}>
    <div onClick={hasChildren ? handleToggle : null} className={classes.item}>
      {icon}
      <Typography className={classes.name}>{name}</Typography>
      {path && <div className={classes.actions}>
        <Checkbox onClick={handleSelect} checked={filter?.has(path) || false} />
      </div>}
    </div>
    <Collapse in={open} className={classes.children}>
      {open && def.quantities.filter(filterProperties).map(def => (
        <Definition key={def.name} def={def} path={childPathPrefix + def.name} />
      ))}
      {open && def.sub_sections.filter(filterProperties).map(def => (
        <Definition key={def.name} name={def.name} path={childPathPrefix + def.name} def={resolveRef(def.sub_section)} />
      ))}
    </Collapse>
  </div>
})
Definition.propTypes = {
  def: PropTypes.object.isRequired,
  name: PropTypes.string,
  path: PropTypes.string
}

const FilterSubMenuArchive = React.memo(({
  value,
  ...rest
}) => {
  return <FilterSubMenu value={value} {...rest}>
    <InputGrid spacing={2}>
      <InputGridItem xs={12}>
        <Definition def={rootSections.find(def => def.name === 'EntryArchive')} />
      </InputGridItem>
    </InputGrid>
  </FilterSubMenu>
})
FilterSubMenuArchive.propTypes = {
  value: PropTypes.string
}

export default FilterSubMenuArchive
