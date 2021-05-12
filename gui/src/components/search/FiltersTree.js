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
import { makeStyles } from '@material-ui/core/styles'
import ArrowForwardIosIcon from '@material-ui/icons/ArrowForwardIos'
import { Typography } from '@material-ui/core'
import FiltersElements from './FiltersElements'

/**
 * Displays the tree-like structure for selecting filters.
 */
const useStyles = makeStyles(theme => ({
  root: {},
  li: {
    display: 'flex',
    alignItems: 'center',
    padding: theme.spacing(0.5),
    borderBottom: `1px solid ${theme.palette.divider}`
  },
  section: {
  },
  link: {
    cursor: 'pointer',
    '&:hover': {
      backgroundColor: theme.palette.secondary
    }
  },
  label: {
    color: theme.palette.text.secondary
  },
  hidden: {
    display: 'none'
  },
  indented: {
    marginLeft: theme.spacing(4)
  },
  arrow: {
    marginLeft: theme.spacing(1),
    fontSize: '0.8rem'
  }
}))

const FiltersTree = React.memo(({
  view,
  onViewChange,
  className
}) => {
  const styles = useStyles()

  // Determine the navigation tree layout
  const tree = useMemo(() => {
    const tree = [
      {
        name: 'Structure',
        children: [
          {name: 'Elements', leaf: true},
          {name: 'Symmetry', leaf: true}
        ]
      },
      {
        name: 'Method',
        children: [
          {
            name: 'Experiment',
            children: [
              {name: 'XPS', leaf: true}
            ]
          },
          {
            name: 'Simulation',
            children: [
              {name: 'DFT', leaf: true},
              {name: 'GW', leaf: true}
            ]
          }
        ]
      },
      {
        name: 'Properties',
        children: [
          {name: 'Electronic', leaf: true},
          {name: 'Vibrational', leaf: true},
          {name: 'Optical', leaf: true}
        ]
      },
      {
        name: 'Metainfo',
        children: [
          {name: 'Origin', leaf: true},
          {name: 'Dataset', leaf: true}
        ]
      }
    ]
    function build(branch, index, level) {
      const name = branch.name
      const children = branch.children
      const leaf = branch.leaf
      const val = <div className={level !== 0 ? styles.indented : undefined} key={index}>
        <div className={clsx(styles.li, leaf && styles.link)} onClick={ leaf ? () => onViewChange(1, name) : () => {} }>
          <Typography className={styles.label} variant="button">{name}</Typography>
          {leaf && <ArrowForwardIosIcon className={styles.arrow}/>}
        </div>
        {children &&
          <div>
            {children.map((child, j) => build(child, j, level + 1))}
          </div>
        }
      </div>
      return val
    }
    return tree.map((branch, i) => build(branch, i, 0))
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  return <div className={clsx(className, styles.root)}>
    <div className={clsx(view !== 'Filters' && styles.hidden)}>{tree}</div>
    <FiltersElements className={clsx(view !== 'Elements' && styles.hidden)}/>
  </div>
})
FiltersTree.propTypes = {
  view: PropTypes.string,
  level: PropTypes.number,
  className: PropTypes.string,
  onViewChange: PropTypes.func
}
FiltersTree.defaultProps = {
  level: 0
}

export default FiltersTree
