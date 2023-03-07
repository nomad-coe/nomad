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
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import { QuantityCell, QuantityRow, QuantityTable } from '../../Quantity'
import { isNil, chunk } from 'lodash'
import {
  Accordion as MuiAccordion,
  AccordionDetails as MuiAccordionDetails,
  AccordionSummary as MuiAccordionSummary,
  Card,
  CardContent,
  CardHeader,
  CardActions,
  Grid,
  Typography,
  makeStyles,
  withStyles
} from '@material-ui/core'
import { traverseDeep } from '../../../utils'

/**
 * Contains components for displaying properties. The main idea is:
 * - Each main property category should be placed inside a PropertyCard.
 * - A property card should have PropertyGrid as a direct descendant.
 * - PropertyItems are used to display individual properties inside a PropertyGrid.
 * - If you wish to subdivide a single PropertyItem into smaller subitems, use a
 *   PropertySubGrid with the set of wanted PropertyItems inside it.
 */

const gridSpacing = 2

/**
 * Card for displaying a set of related properties.
 */
const usePropertyCardStyles = makeStyles(theme => ({
  header: {
    paddingBottom: theme.spacing(1),
    minHeight: 0,
    flexGrow: '0'
  },
  content: {
    height: '100%',
    width: '100%',
    boxSizing: 'border-box',
    minHeight: 0
  },
  action: {
    marginRight: theme.spacing(0),
    marginTop: theme.spacing(0)
  },
  container: {
    display: 'flex',
    flexDirection: 'column',
    height: '100%',
    width: '100%'
  }
}))
export function PropertyCard({className, children, ...headerProps}) {
  const styles = usePropertyCardStyles()

  return <Card className={className} data-testid='property-card'>
    <div className={styles.container}>
      <CardHeader
        {...headerProps}
        className={styles.header}
        classes={{action: styles.action}}
      />
      <CardContent className={styles.content}>
        {children}
      </CardContent>
    </div>
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
const usePropertyCardActionsStyles = makeStyles(theme => ({
  root: {
    marginLeft: 'auto'
  }
}))
export function PropertyCardActions({children}) {
  const styles = usePropertyCardActionsStyles()
  return <CardActions disableSpacing>
    <div className={styles.root}>
      {children}
    </div>
  </CardActions>
}

PropertyCardActions.propTypes = {
  children: PropTypes.any
}

/**
 * For displaying actions at the bottom of a PropertyCard.
 */
const usePropertyTitleStyles = makeStyles(theme => ({
  root: {
    textTransform: 'none',
    fontSize: '0.9rem'
  }
}))
export function PropertyTitle({title, align, className}) {
  const styles = usePropertyTitleStyles()
  return title
    ? <Typography
      variant="button"
      align={align || 'center'}
      className={clsx(className, styles.root)}
    >{title}</Typography>
    : null
}

PropertyTitle.propTypes = {
  title: PropTypes.string,
  align: PropTypes.string,
  className: PropTypes.string
}

/**
 * For displaying a row of properties, typically within a PropertyCard.
 */
export function PropertyGrid({className, children}) {
  return <Grid container spacing={gridSpacing} className={className} style={{height: '100%'}}>
    {children}
  </Grid>
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
    marginBottom: theme.spacing(1)
  },
  column: {
    width: '100%',
    height: '100%',
    display: 'flex',
    flexDirection: 'column'
  },
  content: {
    flex: 1,
    minHeight: 0 // min-height: 0 to allow the item to shrink to fit inside the container.
  }
}))
export function PropertyItem({
  title,
  align,
  className,
  classes,
  children,
  height,
  minHeight,
  ...other
}) {
  const styles = usePropertyItemStyles({classes: classes})
  return <Grid
    item
    {...other}
    className={className}
    style={
      (!isNil(height) || !isNil(minHeight)) && {height: height, minHeight: minHeight}
    }
  >
    <div className={styles.column}>
      <PropertyTitle title={title} className={styles.title}/>
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

/**
 * For displaying the methodology steps for a property.
 */
const usePropertyMethodologyListStyle = makeStyles(theme => ({
  column: {
    width: '100%',
    height: '100%',
    display: 'flex',
    flexDirection: 'column'
  },
  content: {
    flex: 1,
    minHeight: 0 // min-height: 0 to allow the item to shrink to fit inside the container.
  }
}))
export function PropertyMethodologyList({className, classes, children, ...other}) {
  const styles = usePropertyMethodologyListStyle({classes: classes})
  return <Grid item {...other} className={className}>
    <div className={styles.column}>
      <div className={styles.content}>
        {children}
      </div>
    </div>
  </Grid>
}

PropertyMethodologyList.propTypes = {
  className: PropTypes.string,
  classes: PropTypes.object,
  children: PropTypes.any
}

const Accordion = withStyles((theme) => ({
  expanded: {
    margin: theme.spacing(1)
  }
}))(MuiAccordion)

const AccordionSummary = withStyles((theme) => ({
  root: {
    padding: `0 ${theme.spacing(1)}px`,
    minHeight: theme.spacing(4),
    '&$expanded': {
      minHeight: theme.spacing(4)
    }
  },
  content: {
    margin: 0,
    '&$expanded': {
      margin: 0
    }
  },
  expanded: {}
}))(MuiAccordionSummary)

const AccordionDetails = withStyles((theme) => ({
  root: {
    margin: 0,
    padding: `0 ${theme.spacing(1)}px`
  }
}))(MuiAccordionDetails)

/**
 * For displaying methodology steps.
 */
export const PropertyMethodologyItem = React.memo(({title, data, path, columns}) => {
  // Recursively extract items from the data and split the items into equal sized rows
  const rows = useMemo(() => {
    if (!data) return undefined

    const quantities = []
    for (const [key, value] of traverseDeep(data, true)) {
      quantities.push({
        quantity: `${path === '' ? '' : `${path}.`}${key.join('.')}`,
        value: value
      })
    }

    return chunk(quantities, columns)
  }, [data, path, columns])

  if (!data) {
    return null
  }

  return <Accordion defaultExpanded elevation={0}>
    <AccordionSummary IconButtonProps={{size: 'small'}} expandIcon={<ExpandMoreIcon />}>
      <PropertyTitle title={title} />
    </AccordionSummary>
    <AccordionDetails>
      <QuantityTable>
        {rows.map((row) =>
          <QuantityRow key={row[0].quantity}>
            {row.map((column) => <QuantityCell
              key={column.quantity}
              {...column}
            />)}
          </QuantityRow>
        )}
      </QuantityTable>
    </AccordionDetails>
  </Accordion>
})

PropertyMethodologyItem.propTypes = {
  title: PropTypes.string,
  data: PropTypes.object,
  path: PropTypes.string,
  columns: PropTypes.number
}

PropertyMethodologyItem.defaultProps = {
  columns: 4
}
