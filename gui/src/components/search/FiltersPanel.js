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
import clsx from 'clsx'
import { makeStyles } from '@material-ui/core/styles'
import { Paper, Breadcrumbs, Typography } from '@material-ui/core'
import ReloadIcon from '@material-ui/icons/Cached'
import Actions from '../Actions'
import FiltersTree from './FiltersTree'
import FiltersActive from './FiltersActive'

/**
 * Displays the filters panel.
 */
const useStyles = makeStyles(theme => ({
  root: {
    padding: `${theme.spacing(1.5)}px ${theme.spacing(3)}px`,
    boxSizing: 'border-box',
    borderRadius: 0,
    position: 'absolute',
    zIndex: 2,
    display: 'flex',
    flexDirection: 'column',
    minWidth: '100%',
    height: '100%'
  },
  spacer: {
    flexGrow: 1
  },
  breadCrumbs: {
    marginBottom: theme.spacing(2)
  },
  crumb: {
    minHeight: 100,
    cursor: 'pointer',
    color: theme.palette.primary.main
  }
}))

const FiltersPanel = React.memo(({
  className
}) => {
  const styles = useStyles()
  const [breadCrumbs, setBreadCrumbs] = useState(['Filters'])
  const [view, setView] = useState('Filters')

  // List of actionable buttons for the filters
  const actions = [
    {tooltip: 'Clear filters', onClick: () => {}, content: <ReloadIcon/>}
  ]

  // Handling view changes in the tree
  const handleViewChange = useCallback((level, name) => {
    setView(name)
    setBreadCrumbs(old => {
      const nLevels = old.length
      if (level > nLevels - 1) {
        old.push(name)
      } else {
        old = old.slice(0, level + 1)
      }
      return [...old]
    })
  }, [])

  return <Paper elevation={3} className={clsx(className, styles.root)}>
    <Breadcrumbs className={styles.breadCrumbs}>
      { breadCrumbs.map((name, index) => <Typography
        className={styles.crumb}
        variant="button"
        key={index}
        onClick={() => handleViewChange(index, name)}
      >{name}</Typography>) }
    </Breadcrumbs>
    <FiltersTree view={view} onViewChange={handleViewChange}/>
    <FiltersActive/>
    <div className={styles.spacer}/>
    <Actions
      color="primary"
      variant="contained"
      actions={actions}
    />
  </Paper>
})
FiltersPanel.propTypes = {
  className: PropTypes.string
}

export default FiltersPanel
