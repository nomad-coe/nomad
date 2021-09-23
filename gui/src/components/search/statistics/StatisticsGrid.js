
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
import { Paper, Grid } from '@material-ui/core'
import { makeStyles } from '@material-ui/core/styles'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { filterData, useAnchorsValue } from '../SearchContext'

/**
 * Displays a summary of the anchored filter statistics in a grid.
 */
const useStyles = makeStyles(theme => {
  return {
    root: {
      padding: theme.spacing(2),
      paddingTop: theme.spacing(1)
    }
  }
})
const StatisticsGrid = React.memo(({
  className,
  classes
}) => {
  const styles = useStyles(classes)
  const anchors = useAnchorsValue()

  return anchors.length > 0
    ? <Paper className={clsx(className, styles.root)}>
      <Grid container spacing={2}>
        {anchors.map((filter) => {
          const config = filterData[filter].statConfig
          return <Grid item key={filter} {...config.layout}>
            <config.component quantity={filter} visible initialAggSize={10}/>
          </Grid>
        })}
      </Grid>
    </Paper>
    : null
})

StatisticsGrid.propTypes = {
  className: PropTypes.string,
  classes: PropTypes.object
}

export default StatisticsGrid
