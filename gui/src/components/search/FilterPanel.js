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
import React, { useState, useMemo } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { makeStyles, fade } from '@material-ui/core/styles'
import {
  Paper,
  Typography,
  ClickAwayListener
} from '@material-ui/core'
import ClearIcon from '@material-ui/icons/Clear'
import CodeIcon from '@material-ui/icons/Code'
import ArrowBackIcon from '@material-ui/icons/ArrowBack'
import FilterTree from './FilterTree'
import Scrollable from '../visualization/Scrollable'
import Actions from '../Actions'

import FilterElements, { labelElements } from './FilterElements'
import FilterElectronic, { labelElectronic } from './FilterElectronic'
import FilterSymmetry, { labelSymmetry } from './FilterSymmetry'
import FilterDFT, { labelDFT } from './FilterDFT'
import FilterGW, { labelGW } from './FilterGW'
import FilterMethod, { labelMethod } from './FilterMethod'
import FilterMaterial, { labelMaterial } from './FilterMaterial'
import FilterMetadata, { labelMetadata } from './FilterMetadata'
import FilterSummary from './FilterSummary'
import {
  useResetFilters,
  filtersElements,
  filtersMaterial,
  filtersElectronic,
  filtersDFT,
  filtersGW,
  filtersSymmetry,
  filtersMethod,
  filtersMetadata
} from './FilterContext'

/**
 * Displays the filters panel.
 */

const filterTree = [
  {
    name: labelMaterial,
    component: FilterMaterial,
    filters: <FilterSummary quantities={filtersMaterial}/>,
    children: [
      {
        name: labelElements,
        filters: <FilterSummary quantities={filtersElements}/>,
        component: FilterElements
      },
      {
        name: labelSymmetry,
        filters: <FilterSummary quantities={filtersSymmetry}/>,
        component: FilterSymmetry
      }
    ]
  },
  {
    name: labelMethod,
    component: FilterMethod,
    filters: <FilterSummary quantities={filtersMethod}/>,
    children: [
      {
        name: labelDFT,
        filters: <FilterSummary quantities={filtersDFT}/>,
        component: FilterDFT
      },
      {
        name: labelGW,
        filters: <FilterSummary quantities={filtersGW}/>,
        component: FilterGW
      }
    ]
  },
  {
    name: 'Properties',
    children: [
      {
        name: labelElectronic,
        filters: <FilterSummary quantities={filtersElectronic}/>,
        component: FilterElectronic
      }
    ]
  },
  {
    name: labelMetadata,
    filters: <FilterSummary quantities={filtersMetadata}/>,
    component: FilterMetadata
  }
]

const useStyles = makeStyles(theme => {
  // The secondary menu widths are hardcoded. Another option would be to use
  // useLayoutEffect to determine the sizes of the components dynamically, but
  // this seems to be quite a bit less responsive compared to hardcoding the
  // values.
  const padding = theme.spacing(2)
  const widthMedium = 25
  const widthLarge = 42
  return {
    root: {
      boxSizing: 'border-box',
      display: 'flex',
      flexDirection: 'column',
      width: '100%',
      height: '100%'
    },
    header: {
      paddingTop: theme.spacing(0.5),
      paddingBottom: theme.spacing(1),
      paddingLeft: padding,
      paddingRight: padding
    },
    headerSecondary: {
      paddingRight: theme.spacing(0),
      paddingLeft: theme.spacing(0)
    },
    headerText: {
      display: 'flex',
      alignItems: 'center'
    },
    toggles: {
      paddingLeft: padding,
      paddingRight: padding,
      marginBottom: theme.spacing(1),
      height: '2rem'
    },
    toggle: {
      color: fade(theme.palette.action.active, 0.87)
    },
    menuPrimary: {
      zIndex: 3,
      position: 'absolute',
      right: 0,
      top: 0,
      bottom: 0,
      left: 0,
      backgroundColor: theme.palette.background.paper,
      boxSizing: 'border-box'
    },
    menuPrimaryBorder: {
      boxShadow: `1px 0px 0px 0px ${theme.palette.action.selected}`
    },
    // The menu animation uses a transition on the 'transform' property. Notice
    // that animating 'transform' instead of e.g. the 'left' property is much
    // more performant. We also hint the browser that the transform property
    // will be animated using the 'will-change' property: this will pre-optimize
    // the element for animation when possible (the recommendation is to
    // remove/add it when needed, but in this case we keep it on constantly).
    container: {
      zIndex: 2,
      display: 'flex',
      flexDirection: 'column',
      position: 'absolute',
      right: 0,
      top: 0,
      bottom: 0,
      width: `${widthLarge}rem`,
      backgroundColor: theme.palette.background.paper,
      '-webkit-transform': 'none',
      transform: 'none',
      transition: 'transform 250ms',
      flexGrow: 1,
      boxSizing: 'border-box',
      willChange: 'transform'
    },
    containerVisibleLarge: {
      '-webkit-transform': `translateX(${widthLarge}rem)`,
      transform: `translateX(${widthLarge}rem)`
    },
    containerVisibleMedium: {
      '-webkit-transform': `translateX(${widthMedium}rem)`,
      transform: `translateX(${widthMedium}rem)`
    },
    menuSecondary: {
      position: 'absolute',
      right: 0,
      top: 0,
      bottom: 0,
      boxSizing: 'border-box'
    },
    paddingPrimary: {
      paddingTop: `${theme.spacing(1.5)}px`,
      paddingBottom: `${theme.spacing(1.5)}px`
    },
    paddingSecondary: {
      padding: `${theme.spacing(1.5)}px ${padding}px`
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
  isMenuOpen,
  resultType,
  className,
  onResultTypeChange,
  onIsMenuOpenChange
}) => {
  const styles = useStyles()
  const [view, setView] = useState('Filters')
  const resetFilters = useResetFilters()

  const views = useMemo(() => {
    const viewList = []
    function addView(item) {
      // We keep all views active and keep rendering them even if they are not
      // visible. This way the application feels more responsive when switching
      // between panels, with the cost of the application being less responsive
      // on changing individual values due to renders on hidden components.
      if (item.component) {
        viewList.push(<item.component
          key={viewList.length}
          visible={view === item.name}
          className={clsx(view !== item.name && styles.menuHidden)}
        ></item.component>)
      }
      let children = item.children
      if (children) {
        for (let child of children) {
          addView(child)
        }
      }
    }
    for (let item of filterTree) {
      addView(item)
    }
    return viewList
  }, [styles, view])

  // Primary menu actions
  const actionsPrimary = useMemo(() => (
    [
    // {
    //   tooltip: 'View the API call for the selected filters',
    //   content: <CodeIcon/>,
    //   onClick: () => {}
    // },
      {
        tooltip: 'Clear filters',
        content: <ClearIcon/>,
        onClick: () => resetFilters()
      }
    ]
  // eslint-disable-next-line react-hooks/exhaustive-deps
  ), [])

  // Secondary menu actions
  const actionsSecondary = useMemo(() => (
    [{
      tooltip: 'Hide filter panel',
      content: <ArrowBackIcon/>,
      onClick: () => { onIsMenuOpenChange(false) }
    }]
  // eslint-disable-next-line react-hooks/exhaustive-deps
  ), [])

  // return <ClickAwayListener onClickAway={() => onIsMenuOpenChange(false)}>
  return <ClickAwayListener onClickAway={() => {}}>
    <div className={clsx(className, styles.root)}>
      <Scrollable className={clsx(styles.menuPrimary, isMenuOpen && styles.menuPrimaryBorder)}>
        <div className={styles.paddingPrimary}>
          <Actions
            header={<Typography className={styles.headerText} variant="button">Filters</Typography>}
            variant="icon"
            actions={actionsPrimary}
            className={styles.header}
          />
          <FilterTree
            filterTree={filterTree}
            view={view}
            isMenuOpen={isMenuOpen}
            resultType={resultType}
            onResultTypeChange={onResultTypeChange}
            onViewChange={setView}
            onIsMenuOpenChange={onIsMenuOpenChange}
          />
        </div>
      </Scrollable>
      <Paper
        elevation={4}
        className={clsx(styles.container, isMenuOpen && (view !== labelElements ? styles.containerVisibleMedium : styles.containerVisibleLarge))}
      >
        <div className={clsx(styles.menuSecondary, view !== labelElements ? styles.menuMedium : styles.menuLarge)}>
          <Scrollable>
            <div className={styles.paddingSecondary}>
              <Actions
                header={<Typography className={styles.headerText} variant="button">{view}</Typography>}
                variant="icon"
                actions={actionsSecondary}
                className={clsx(styles.header, styles.headerSecondary)}
              />
              {views}
            </div>
          </Scrollable>
        </div>
      </Paper>
    </div>
  </ClickAwayListener>
})
FilterPanel.propTypes = {
  isMenuOpen: PropTypes.bool,
  resultType: PropTypes.string.isRequired,
  className: PropTypes.string,
  onResultTypeChange: PropTypes.func,
  onIsMenuOpenChange: PropTypes.func
}

export default FilterPanel
