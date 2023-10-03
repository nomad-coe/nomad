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
import React, { useCallback, useEffect, useMemo, useState, useContext } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { FilterSubMenu, filterMenuContext } from './FilterMenu'
import { InputGrid, InputGridItem } from '../input/InputGrid'
import { Box, Collapse, makeStyles } from '@material-ui/core'
import { SectionMDef, SubSectionMDef, useGlobalMetainfo } from '../../archive/metainfo'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ChevronRightIcon from '@material-ui/icons/ChevronRight'
import { useSearchContext } from '../SearchContext'
import InputItem from '../input/InputItem'
import { InputTextQuantity } from '../input/InputText'
import { useErrors } from '../../errors'
import { useDataStore } from '../../DataStore'
import InputField from '../input/InputField'
import InputRadio from '../input/InputRadio'
import { useApi } from '../../api'
import InputHeader from '../input/InputHeader'

const filterProperties = def => !(def.name.startsWith('x_') || def.virtual)

const useDefinitionStyles = makeStyles(theme => ({
  root: {},
  item: {
    '&:hover': {
      backgroundColor: theme.palette.grey[100]
    },
    cursor: 'pointer',
    width: '100%',
    display: 'flex',
    flexDirection: 'row',
    alignItems: 'center',
    flexWrap: 'nowrap',
    paddingRight: theme.spacing(1),
    boxSizing: 'border-box'
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

const Definition = React.memo(function Definition({node, name, path, className}) {
  const classes = useDefinitionStyles()
  const [open, setOpen] = useState(!path)
  const { useFilterState } = useSearchContext()
  const [filter, setFilter] = useFilterState('quantities')

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

  const children = Object.keys(node)
  const hasChildren = children.length > 0

  const icon = useMemo(() => {
    const iconProps = {
      className: classes.icon,
      fontSize: 'small'
    }
    return hasChildren
      ? (open ? <ExpandMoreIcon {...iconProps} /> : <ChevronRightIcon {...iconProps}/>)
      : <div className={classes.icon} />
  }, [hasChildren, open, classes.icon])

  const childPathPrefix = path + '.'
  const indeterminate = useMemo(() => {
    if (!filter) {
      return false
    }
    return hasChildren && Array.from(filter).find(fullPath => fullPath.startsWith(childPathPrefix)) !== undefined
  }, [filter, hasChildren, childPathPrefix])

  return <div className={clsx(className, classes.root)}>
    <div onClick={hasChildren ? handleToggle : null} className={classes.item}>
      {icon}
      <InputItem
        value={path}
        label={name}
        selected={filter?.has(path) || false}
        onChange={handleSelect}
        variant="checkbox"
        labelPlacement="start"
        disableStatistics
        disableLabelClick
        disableSelect={!path}
        indeterminate={indeterminate}
      />
    </div>
    {hasChildren && (
      <Collapse in={open} className={classes.children}>
        {open && children.map((name, index) => (
          <Definition
            key={index}
            name={name}
            node={node[name]}
            path={childPathPrefix + name} />
        ))}
      </Collapse>
    )}
  </div>
})
Definition.propTypes = {
  node: PropTypes.object.isRequired,
  name: PropTypes.string,
  path: PropTypes.string,
  className: PropTypes.string
}

const FilterSubMenuMetadata = React.memo(({
  id,
  ...rest
}) => {
  const {api} = useApi()
  const {selected, open} = useContext(filterMenuContext)
  const visible = open && id === selected
  const authenticated = api?.keycloak?.authenticated
  const dataStore = useDataStore()
  const globalMetainfo = useGlobalMetainfo()
  const [[options, tree], setOptions] = useState([[], {}])
  const {raiseError} = useErrors()
  useEffect(() => {
    if (!globalMetainfo || !visible) {
      return
    }
    const root = globalMetainfo.getEntryArchiveDefinition()
    const subSections = root?.sub_sections
    const tree = {}
    if (!subSections) {
      return []
    }
    const options = []
    const defsSet = new Set()
    function addDef(def, prefix, node) {
      const fullName = prefix ? `${prefix}.${def.name}` : def.name
      node[def.name] = node[def.name] || {}
      node = node[def.name]
      if (defsSet.has(def)) {
        return
      }
      defsSet.add(def)
      options.push({value: fullName})
      if (def.m_def === SubSectionMDef && def.sub_section) {
        def = def.sub_section
      }
      if (def.m_def === SectionMDef) {
        def._allProperties.filter(filterProperties).forEach(def => addDef(def, fullName, node))
        dataStore.getAllInheritingSections(def).forEach(def => {
          def._allProperties.filter(filterProperties).forEach(def => addDef(def, fullName, node))
        })
      }
    }
    subSections.forEach(def => addDef(def, null, tree))
    setOptions([options, tree])
  }, [raiseError, dataStore, globalMetainfo, setOptions, visible])

  return <FilterSubMenu id={id} {...rest}>
    <InputGrid>
      <InputGridItem xs={12}>
        <InputRadio
          quantity="visibility"
          label="Visibility"
          initialValue={authenticated ? 'visible' : 'public'}
          options={{
            all: {label: 'All', disabled: !authenticated, tooltip: 'Consider all entries.'},
            public: {label: 'Public', disabled: false, tooltip: 'Consider all entries that can be publically downloaded, i.e. only published entries without embargo.'},
            visible: {label: 'Visible', disabled: !authenticated, tooltip: 'Consider all entries that are visible to you. This includes entries with embargo or unpublished entries that belong to you or are shared with you.'},
            shared: {label: 'Shared', disabled: !authenticated, tooltip: 'Only consider entries that belong to you or are shared with you.'},
            user: {label: 'User', disabled: !authenticated, tooltip: 'Only consider entries that belong to you.'},
            staging: {label: 'Unpublished', disabled: !authenticated, tooltip: 'Only search through unpublished entries.'}
          }}
        ></InputRadio>
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputField
          quantity="entry_id"
          visible={visible}
          disableStatistics
          disableOptions
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputField
          quantity="upload_id"
          visible={visible}
          disableStatistics
          disableOptions
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputField
          quantity="upload_name"
          visible={visible}
          disableStatistics
          disableOptions
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputField
          quantity="results.material.material_id"
          visible={visible}
          disableStatistics
          disableOptions
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputField
          quantity="datasets.dataset_id"
          visible={visible}
          disableStatistics
          disableOptions
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputHeader quantity="quantities" disableStatistics/>
        <Box marginBottom={1}>
          <InputTextQuantity
            quantity="quantities"
            suggestions={options}
          />
        </Box>
        {Object.keys(tree).map((name, index) => (
          <Definition node={tree[name]} name={name} path={name} key={index} />
        ))}
      </InputGridItem>
    </InputGrid>
  </FilterSubMenu>
})
FilterSubMenuMetadata.propTypes = {
  id: PropTypes.string
}

export default FilterSubMenuMetadata
