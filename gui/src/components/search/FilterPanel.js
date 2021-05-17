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
  Typography,
  ClickAwayListener
} from '@material-ui/core'
import ClearIcon from '@material-ui/icons/Clear'
import CodeIcon from '@material-ui/icons/Code'
import FiltersTree from './FilterTree'
import FilterElements from './FilterElements'
import FilterSymmetry from './FilterSymmetry'
import Actions from '../Actions'

/**
 * Displays the filters panel.
 */

const useStyles = makeStyles(theme => {
  const padding = theme.spacing(2)
  const widthMedium = 30
  const widthLarge = 42
  return {
    root: {
      boxSizing: 'border-box',
      borderRadius: 0,
      position: 'absolute',
      zIndex: 2,
      display: 'flex',
      flexDirection: 'column',
      minWidth: '100%',
      height: '100%'
    },
    header: {
      paddingTop: theme.spacing(1),
      paddingBottom: theme.spacing(2),
      paddingLeft: padding,
      paddingRight: padding
    },
    menuPrimary: {
      zIndex: 1,
      backgroundColor: theme.palette.background.paper,
      paddingTop: theme.spacing(1.5),
      paddingBottom: theme.spacing(1.5),
      flexGrow: 1,
      borderRight: `1px solid ${theme.palette.divider}`
    },
    // The secondary menu widths are hardcoded. Another option would be to use
    // useLayoutEffect to determine the sizes of the components dynamically, but
    // this seems to be quite a bit less responsive compared to hardcoding the
    // values.
    container: {
      position: 'absolute',
      zIndex: 0,
      right: 0,
      left: 0,
      top: 0,
      bottom: 0,
      transition: 'right 200ms',
      flexGrow: 1,
      boxSizing: 'border-box'
    },
    menuSecondary: {
      position: 'absolute',
      right: 0,
      top: 0,
      bottom: 0,
      padding: `${theme.spacing(1.5)}px ${padding}px`,
      boxSizing: 'border-box'
    },
    containerLarge: {
      right: `-${widthLarge}rem`
    },
    containerMedium: {
      right: `-${widthMedium}rem`
    },
    menuLarge: {
      width: `${widthLarge}rem`
    },
    menuMedium: {
      width: `${widthMedium}rem`
    },
    menuHidden: {
      display: 'none'
    }
  }
})

const FilterPanel = React.memo(({
  resultType,
  className,
  onResultTypeChange
}) => {
  const styles = useStyles()
  const [view, setView] = useState('Filters')
  const [isMenuVisible, setIsMenuVisible] = useState(false)

  // Handling view changes in the tree
  const handleViewChange = useCallback((level, name) => {
    setView(name)
    setIsMenuVisible(name !== 'filters')
  }, [])

  const actions = useMemo(() => (
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
  ), [handleViewChange])

  return <ClickAwayListener onClickAway={() => setIsMenuVisible(false)}>
    <div className={clsx(className, styles.root)}>
      <div className={styles.menuPrimary}>
        <Actions
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
      </div>
      <Paper
        elevation={4}
        className={clsx(styles.container, isMenuVisible && (view !== 'Elements / Formula' ? styles.containerMedium : styles.containerLarge))}
      >
        <div
          className={clsx(styles.menuSecondary, isMenuVisible && (view !== 'Elements / Formula' ? styles.menuMedium : styles.menuLarge))}
        >
          <FilterElements className={clsx(styles.menu, view !== 'Elements / Formula' && styles.menuHidden)}/>
          <FilterSymmetry className={clsx(styles.menu, view !== 'Symmetry / Prototypes' && styles.menuHidden)}/>
        </div>
      </Paper>
    </div>
  </ClickAwayListener>
})
FilterPanel.propTypes = {
  resultType: PropTypes.string.isRequired,
  className: PropTypes.string,
  onResultTypeChange: PropTypes.func
}

export default FilterPanel
