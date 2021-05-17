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
import React, { useCallback, useState, useMemo } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { makeStyles } from '@material-ui/core/styles'
import {
  Paper,
  Breadcrumbs,
  Typography
} from '@material-ui/core'
import ArrowBackIcon from '@material-ui/icons/ArrowBack'
import ClearIcon from '@material-ui/icons/Clear'
import CodeIcon from '@material-ui/icons/Code'
import FiltersTree from './FiltersTree'
import Actions from '../Actions'

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
  header: {
    paddingTop: theme.spacing(1),
    paddingBottom: theme.spacing(2)
  },
  breadCrumbs: {
    marginBottom: theme.spacing(1.5)
  },
  back: {
    display: 'flex',
    alignItems: 'center'
  },
  crumb: {
    cursor: 'pointer',
    color: theme.palette.primary.main
  }
}))

const FiltersPanel = React.memo(({
  resultType,
  className,
  onResultTypeChange
}) => {
  const styles = useStyles()
  const [breadCrumbs, setBreadCrumbs] = useState(['Filters'])
  const [view, setView] = useState('Filters')

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

  const actions = useMemo(() => (
    // [{
    //   tooltip: 'Go back to main menu',
    //   onClick: () => { handleViewChange(0, 'Filters') },
    //   content: <div className={styles.back}><ArrowBackIcon/>BACK</div>
    // },
    [{
      tooltip: 'View the API call for the selected filters',
      content: <CodeIcon/>,
      onClick: () => { handleViewChange(0, 'Filters') }
    },
    {
      tooltip: 'Clear filters',
      content: <ClearIcon/>,
      onClick: () => { handleViewChange(0, 'Filters') }
    }]
  ), [handleViewChange, styles, view])

  console.log(resultType)
  return <Paper elevation={3} className={clsx(className, styles.root)}>
    <Actions
      // header={<Breadcrumbs>
      //   { breadCrumbs.map((name, index) => <Typography
      //     className={styles.crumb}
      //     variant="button"
      //     key={index}
      //     onClick={() => handleViewChange(index, name)}
      //   >{name}</Typography>) }
      // </Breadcrumbs>
      // }
      header={<Typography
        variant="button"
      >Filters
      </Typography>}
      variant="icon"
      actions={actions}
      className={styles.header}
    />
    <FiltersTree
      view={view}
      resultType={resultType}
      onResultTypeChange={onResultTypeChange}
      onViewChange={handleViewChange}
    />
    <div className={styles.spacer}/>
  </Paper>
})
FiltersPanel.propTypes = {
  resultType: PropTypes.string.isRequired,
  className: PropTypes.string,
  onResultTypeChange: PropTypes.func
}

export default FiltersPanel
