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
import React, { useState, useCallback, useContext, useEffect } from 'react'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import {
  List,
  ListItem,
  ListItemText,
  ListItemIcon,
  Divider,
  Paper,
  Typography
} from '@material-ui/core'
import ArrowBackIcon from '@material-ui/icons/ArrowBack'
import NavigateNextIcon from '@material-ui/icons/NavigateNext'
import ClearIcon from '@material-ui/icons/Clear'
import Scrollable from '../../visualization/Scrollable'
import FilterSummary from '../FilterSummary'
import { Actions, Action } from '../../Actions'
import { quantityGroups, useResetFilters } from '../FilterContext'

const useFilterMenuStyles = makeStyles(theme => {
  const padding = theme.spacing(2)
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
      paddingBottom: theme.spacing(0.5),
      paddingLeft: padding,
      paddingRight: padding
    },
    headerText: {
      display: 'flex',
      alignItems: 'center'
    },
    menu: {
      zIndex: 3,
      position: 'absolute',
      right: 0,
      top: 0,
      bottom: 0,
      left: 0,
      backgroundColor: theme.palette.background.paper,
      boxSizing: 'border-box'
    },
    menuBorder: {
      boxShadow: `1px 0px 0px 0px ${theme.palette.action.selected}`
    },
    padding: {
      paddingTop: `${theme.spacing(1.5)}px`,
      paddingBottom: `${theme.spacing(1.5)}px`
    },
    button: {
      marginRight: 0
    }
  }
})

export const FilterMenu = React.memo(({
  className,
  children
}) => {
  const styles = useFilterMenuStyles()
  const { open } = useContext(filterMenuContext)
  const resetFilters = useResetFilters()

  // Unfortunately the ClickAwayListener does not play nicely together with
  // Menus/Select/Popper. When using Portals, the clicks are registered wrong.
  // When Portals are disabled (disablePortal), their positioning goes haywire.
  // The clicks outside are thus detected by individual event listeners that
  // toggle the menu state.
  return <div className={clsx(className, styles.root)}>
    <Scrollable className={clsx(styles.menu, open && styles.menuBorder)}>
      <div className={styles.padding}>
        <Actions
          header={<Typography className={styles.headerText} variant="button">Filters</Typography>}
          className={styles.header}
        >
          <Action
            tooltip="Clear filters"
            onClick={() => resetFilters()}
            className={styles.button}
          >
            <ClearIcon/>
          </Action>
        </Actions>
        <List dense className={styles.list}>
          <Divider/>
          {children}
        </List>
      </div>
    </Scrollable>
  </div>
})
FilterMenu.propTypes = {
  className: PropTypes.string,
  children: PropTypes.node
}

const useFilterMenuItemStyles = makeStyles(theme => {
  return {
    root: {},
    label: {
      textTransform: 'capitalize'
    },
    listIcon: {
      fontsize: '1rem',
      minWidth: '1.5rem'
    },
    arrow: {
      marginLeft: theme.spacing(1),
      fontSize: '1.5rem'
    },
    gutters: {
      paddingLeft: theme.spacing(3.0),
      paddingRight: theme.spacing(2.35)
    },
    listItem: {
      position: 'relative',
      height: '2.6rem'
    },
    divider: {
      width: '100%'
    }
  }
})
export const FilterMenuItem = React.memo(({
  value,
  onClick,
  disableButton,
  depth
}) => {
  const styles = useFilterMenuItemStyles()
  const theme = useTheme()
  const group = quantityGroups.get(value)
  const { selected, open, onChange } = useContext(filterMenuContext)
  const handleClick = disableButton ? undefined : (onClick || onChange)
  const opened = open && value === selected

  return <>
    <ListItem
      button={!!handleClick}
      className={styles.listItem}
      classes={{gutters: styles.gutters}}
      onClick={handleClick && (() => handleClick(value))}
    >
      <ListItemText
        style={{marginLeft: theme.spacing(depth * 1.8)}}
        primaryTypographyProps={{
          color: opened ? 'primary' : 'initial',
          className: styles.label
        }}
        primary={value}
      />
      {handleClick && <ListItemIcon className={styles.listIcon}>
        <NavigateNextIcon color={opened ? 'primary' : 'action'} className={styles.arrow}/>
      </ListItemIcon>}
    </ListItem>
    {group && <FilterSummary quantities={group}/>}
    <Divider className={styles.divider}/>
  </>
})

FilterMenuItem.propTypes = {
  value: PropTypes.string,
  onClick: PropTypes.func,
  disableButton: PropTypes.bool,
  depth: PropTypes.number
}
FilterMenuItem.defaultProps = {
  depth: 0
}

export const filterMenuContext = React.createContext()
export const FilterMenuContext = React.memo(({
  selected,
  onSelectedChange,
  open,
  onOpenChange,
  children
}) => {
  const [size, setSize] = useState('medium')

  const handleChange = useCallback((newValue) => {
    if (newValue !== selected) {
      onOpenChange(true)
    } else {
      onOpenChange(old => !old)
    }
    onSelectedChange && onSelectedChange(newValue)
  }, [selected, onSelectedChange, onOpenChange])

  return <filterMenuContext.Provider value={{
    selected: selected,
    open: open,
    onOpenChange: onOpenChange,
    size: size,
    onSizeChange: setSize,
    onChange: handleChange
  }}>
    {children}
  </filterMenuContext.Provider>
})

FilterMenuContext.propTypes = {
  selected: PropTypes.string,
  onSelectedChange: PropTypes.func,
  open: PropTypes.bool,
  onOpenChange: PropTypes.func,
  children: PropTypes.node
}

const useFilterSubMenuStyles = makeStyles(theme => ({
  root: {
    width: '100%'
  },
  hidden: {
    display: 'none'
  }
}))

const useFilterSubMenusStyles = makeStyles(theme => {
  // The secondary menu widths are hardcoded. Another option would be to use
  // useLayoutEffect to determine the sizes of the components dynamically, but
  // this seems to be quite a bit less responsive compared to hardcoding the
  // values.
  const padding = theme.spacing(2)
  const widthMedium = 25
  const widthLarge = 43
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
      paddingBottom: theme.spacing(1.5),
      paddingRight: theme.spacing(0),
      paddingLeft: theme.spacing(0)
    },
    headerText: {
      display: 'flex',
      alignItems: 'center'
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
    containerMedium: {
      '-webkit-transform': `translateX(${widthMedium}rem)`,
      transform: `translateX(${widthMedium}rem)`
    },
    containerLarge: {
      '-webkit-transform': `translateX(${widthLarge}rem)`,
      transform: `translateX(${widthLarge}rem)`
    },
    menu: {
      position: 'absolute',
      right: 0,
      top: 0,
      bottom: 0,
      boxSizing: 'border-box'
    },
    padding: {
      padding: `${theme.spacing(1.5)}px ${padding}px`
    },
    menuMedium: {
      width: `${widthMedium}rem`
    },
    menuLarge: {
      width: `${widthLarge}rem`
    },
    button: {
      marginRight: 0
    }
  }
})

export const FilterSubMenus = React.memo(({
  children
}) => {
  const styles = useFilterSubMenusStyles()
  const { selected, open, onOpenChange, size } = useContext(filterMenuContext)
  const [menuStyle, containerStyle] = {
    medium: [styles.menuMedium, styles.containerMedium],
    large: [styles.menuLarge, styles.containerLarge]
  }[size]

  return <Paper
    elevation={4}
    className={clsx(styles.container, open && containerStyle)}
  >
    <div className={clsx(styles.menu, menuStyle)}>
      <Scrollable>
        <div className={styles.padding}>
          <Actions
            header={<Typography className={styles.headerText} variant="button">{selected}</Typography>}
            className={styles.header}
          >
            <Action
              tooltip="Hide filter panel"
              onClick={() => { onOpenChange(false) }}
              className={styles.button}
            >
              <ArrowBackIcon/>
            </Action>
          </Actions>
          {children}
        </div>
      </Scrollable>
    </div>
  </Paper>
})
FilterSubMenus.propTypes = {
  sizes: PropTypes.object,
  children: PropTypes.node
}

export const FilterSubMenu = React.memo(({
  value,
  size,
  children
}) => {
  const styles = useFilterSubMenuStyles()
  const { selected, onSizeChange } = useContext(filterMenuContext)
  const visible = value === selected
  useEffect(() => {
    if (visible) {
      onSizeChange(size)
    }
  }, [size, visible, onSizeChange])

  return <div className={clsx(styles.root, !visible && styles.hidden)}>
    {children}
  </div>
})
FilterSubMenu.propTypes = {
  value: PropTypes.string,
  size: PropTypes.oneOf(['medium', 'large']),
  children: PropTypes.node
}
FilterSubMenu.defaultProps = {
  size: 'medium'
}

export default FilterSubMenu
