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
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { makeStyles, fade } from '@material-ui/core/styles'
import {
  List,
  ListItem,
  ListItemText,
  ListItemIcon,
  ListSubheader,
  Divider
} from '@material-ui/core'
import {
  ToggleButton,
  ToggleButtonGroup
} from '@material-ui/lab'
import Checkbox from '@material-ui/core/Checkbox'
import NavigateNextIcon from '@material-ui/icons/NavigateNext'
import FiltersElements from './FiltersElements'

/**
 * Displays the tree-like structure for selecting filters.
 */
const useStyles = makeStyles(theme => ({
  root: {},
  list: {
    paddingTop: 0,
    paddingBottom: 0
  },
  section: {
  },
  hidden: {
    display: 'none'
  },
  listHeader: {
    paddingLeft: 0,
    height: '2.5rem',
    lineHeight: '2.5rem',
    color: theme.palette.text.primary
  },
  listIcon: {
    minWidth: '1.5rem'
  },
  listItem: {
    height: '2.5rem'
  },
  arrow: {
    marginLeft: theme.spacing(1),
    fontSize: '1.5rem'
  },
  toggles: {
    marginBottom: theme.spacing(1),
    height: '2rem'
  },
  toggle: {
    color: fade(theme.palette.action.active, 0.87)
  },
  gutters: {
    paddingLeft: '0.5rem',
    paddingRight: '0.5rem'
  },
  selected: {
    '&$selected': {
      backgroundColor: theme.palette.secondary.main,
      color: 'white'
    }
  }
}))

const FiltersTree = React.memo(({
  resultType,
  view,
  onViewChange,
  onResultTypeChange,
  className
}) => {
  const styles = useStyles()

  // Determine the navigation tree layout
  const tree = useMemo(() => {
    const tree = [
      {
        name: 'Structure',
        children: [
          {name: 'Elements'},
          {name: 'Classification'},
          {name: 'Symmetry'}
        ]
      },
      {
        name: 'Method',
        children: [
          {name: 'DFT', onChecked: () => {}},
          {name: 'GW', onChecked: () => {}},
          {name: 'XPS', onChecked: () => {}}
        ]
      },
      {
        name: 'Properties',
        children: [
          {name: 'Electronic'},
          {name: 'Vibrational'},
          {name: 'Optical'}
        ]
      },
      {
        name: 'Metainfo',
        children: [
          {name: 'Origin'},
          {name: 'Dataset'}
        ]
      }
    ]
    function buildList(branch, i) {
      const name = branch.name
      const children = branch.children
      return <List
        key={i}
        dense
        className={styles.list}
        subheader={
          <ListSubheader
            component="div"
            id="nested-list-subheader"
            className={styles.listHeader}
          >
            {name}
          </ListSubheader>
        }
      >
        <Divider/>
        {children.map((child, j) => {
          const childŃame = child.name
          return <ListItem
            divider
            button
            key={j}
            className={styles.listItem}
            classes={{gutters: styles.gutters}}
          >
            <ListItemIcon
              className={styles.listIcon}
            >
              {child.onChecked &&
                <Checkbox
                  edge="start"
                  size="small"
                  checked={true}
                  onChange={(event) => {
                    console.log('Check')
                  }}
                  tabIndex={-1}
                  disableRipple
                />}
            </ListItemIcon>
            <ListItemText
              primary={childŃame}
              onClick={() => onViewChange(1, childŃame)}
            />
            <ListItemIcon
              className={styles.listIcon}
              onClick={() => onViewChange(1, childŃame)}
            >
              <NavigateNextIcon className={styles.arrow}/>
            </ListItemIcon>
          </ListItem>
        })}
      </List>
    }
    return tree.map((section, index) => buildList(section, index))
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  return <div className={clsx(className, styles.root)}>
    <div
      className={clsx(view !== 'Filters' && styles.hidden)}
    >
      <ToggleButtonGroup
        size="small"
        exclusive
        value={resultType}
        onChange={onResultTypeChange}
      >
        <ToggleButton
          value="entries"
          classes={{root: styles.toggle, selected: styles.selected}}
        >Entries
        </ToggleButton>
        <ToggleButton
          value="materials"
          classes={{root: styles.toggle, selected: styles.selected}}
        >Materials
        </ToggleButton>
        <ToggleButton
          value="datasets"
          classes={{root: styles.toggle, selected: styles.selected}}
        >Datasets
        </ToggleButton>
        <ToggleButton
          value="uploads"
          classes={{root: styles.toggle, selected: styles.selected}}
        >Uploads
        </ToggleButton>
      </ToggleButtonGroup>
      {tree}
    </div>
    <FiltersElements className={clsx(view !== 'Elements' && styles.hidden)}/>
  </div>
})
FiltersTree.propTypes = {
  resultType: PropTypes.string.isRequired,
  view: PropTypes.string,
  level: PropTypes.number,
  className: PropTypes.string,
  onViewChange: PropTypes.func,
  onResultTypeChange: PropTypes.func
}
FiltersTree.defaultProps = {
  level: 0
}

export default FiltersTree
