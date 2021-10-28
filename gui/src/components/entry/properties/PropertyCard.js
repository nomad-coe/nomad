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
import {
  Card,
  CardContent,
  CardHeader,
  Grid,
  Typography,
  makeStyles
} from '@material-ui/core'

/**
 * Card for displaying a set of related properties.
 */
export function PropertyCard({children, ...headerProps}) {
  return <Card>
    <CardHeader {...headerProps} />
    <CardContent>
      {children}
    </CardContent>
  </Card>
}
PropertyCard.propTypes = {
  children: PropTypes.oneOfType([
    PropTypes.arrayOf(PropTypes.node),
    PropTypes.node
  ])
}

/**
 * For displaying a row of properties, typically within a PropertyCard.
 */
export function PropertyGrid({children}) {
  return <Grid container spacing={2}>
    {children}
  </Grid>
}

PropertyGrid.propTypes = {
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
export function PropertyItem({title, classes, children, height, ...other}) {
  const styles = usePropertyItemStyles({classes: classes})
  return <Grid item {...other} style={height && {height: height}}>
    <div className={styles.column}>
      {title && <Typography variant="button" align='center' className={styles.title}>{title}</Typography>}
      <div className={styles.content}>
        {children}
      </div>
    </div>
  </Grid>
}

PropertyItem.propTypes = {
  title: PropTypes.string,
  height: PropTypes.string.isRequired,
  classes: PropTypes.object,
  children: PropTypes.any
}

PropertyItem.defaultProps = {
  height: '400px'
}
