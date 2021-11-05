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
import React from 'react'
import PropTypes from 'prop-types'
import { isNil } from 'lodash'
import {
  Card,
  CardContent,
  CardHeader,
  CardActions,
  Grid,
  Typography,
  makeStyles
} from '@material-ui/core'

const gridSpacing = 2

/**
 * Card for displaying a set of related properties.
 */
const usePropertyCardStyles = makeStyles(theme => ({
  header: {
    paddingBottom: theme.spacing(1)
  },
  action: {
    marginRight: theme.spacing(0),
    marginTop: theme.spacing(0)
  }
}))
export function PropertyCard({className, children, ...headerProps}) {
  const styles = usePropertyCardStyles()

  return <Card className={className}>
    <CardHeader {...headerProps} className={styles.header} classes={{action: styles.action}}/>
    {children}
  </Card>
}
PropertyCard.propTypes = {
  className: PropTypes.string,
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
}

/**
 * For displaying actions at the bottom of a PropertyCard.
 */
export function PropertyCardActions({children}) {
  return <CardActions disableSpacing>
    <div style={{marginLeft: 'auto'}}>
      {children}
    </div>
  </CardActions>
}

PropertyCardActions.propTypes = {
  children: PropTypes.any
}

/**
 * For displaying a row of properties, typically within a PropertyCard.
 */
export function PropertyGrid({className, children}) {
  return <CardContent>
    <Grid container spacing={gridSpacing} className={className}>
      {children}
    </Grid>
  </CardContent>
}

PropertyGrid.propTypes = {
  className: PropTypes.string,
  children: PropTypes.any
}

/**
 * For displaying a grid of PropertyItems inside another PropertyItem.
 */
const usePropertySubGridStyles = makeStyles(theme => ({
  root: {
    height: `calc(100% + ${theme.spacing(gridSpacing)}px)`
  }
}))
export function PropertySubGrid({children}) {
  const styles = usePropertySubGridStyles()
  return <Grid container spacing={gridSpacing} className={styles.root}>
    {children}
  </Grid>
}

PropertySubGrid.propTypes = {
  className: PropTypes.string,
  children: PropTypes.any
}

/**
 * For displaying an individual property, typically within a PropertyRow.
 */
const usePropertyItemStyles = makeStyles(theme => ({
  title: {
    marginBottom: theme.spacing(1),
    textTransform: 'none',
    fontSize: '0.9rem'
  },
  column: {
    width: '100%',
    height: '100%',
    display: 'flex',
    flexDirection: 'column'
  },
  content: {
    flex: 1,
    minHeight: 0 // added min-height: 0 to allow the item to shrink to fit inside the container.
  }
}))
export function PropertyItem({title, align, className, classes, children, height, minHeight, ...other}) {
  const styles = usePropertyItemStyles({classes: classes})
  return <Grid item {...other} className={className} style={(!isNil(height) || !isNil(minHeight)) && {height: height, minHeight: minHeight}}>
    <div className={styles.column}>
      {title && <Typography variant="button" align={align || 'center'} className={styles.title}>{title}</Typography>}
      <div className={styles.content}>
        {children}
      </div>
    </div>
  </Grid>
}

PropertyItem.propTypes = {
  title: PropTypes.string,
  align: PropTypes.string,
  height: PropTypes.string.isRequired,
  minHeight: PropTypes.string,
  className: PropTypes.string,
  classes: PropTypes.object,
  children: PropTypes.any
}

PropertyItem.defaultProps = {
  height: '400px'
}
