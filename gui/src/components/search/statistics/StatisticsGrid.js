
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
import { Paper } from '@material-ui/core'
import { makeStyles } from '@material-ui/core/styles'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { MuuriComponent, getResponsiveStyle } from 'muuri-react'
import { filterData, useAnchorsValue, widthMapping } from '../SearchContext'

/**
 * Displays a summary of the anchored filter statistics in a grid.
 */
const useStyles = makeStyles(theme => {
  return {
    root: {
      padding: theme.spacing(1),
      paddingTop: theme.spacing(0)
    },
    muuriInnerItem: {
      backgroundColor: theme.palette.background.paper,
      position: 'absolute',
      left: theme.spacing(1),
      top: theme.spacing(1),
      right: theme.spacing(1),
      bottom: theme.spacing(1)
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
      <MuuriComponent
        dragEnabled
        dragHandle=".dragHandle"
      >
        {anchors.map((filter) => {
          const config = filterData[filter].statConfig
          const layout = config.layout
          const responsiveStyle = getResponsiveStyle({
            columns: widthMapping[layout.widthDefault]['lg'] / 12,
            margin: 0,
            ratio: layout.ratioDefault
          })
          return <div key={filter} style={responsiveStyle}>
            <div className={styles.muuriInnerItem}>
              <config.component
                quantity={filter}
                visible
                initialAggSize={10}
                draggable
              />
            </div>
          </div>
        })}
      </MuuriComponent>
    </Paper>
    : null
})

StatisticsGrid.propTypes = {
  className: PropTypes.string,
  classes: PropTypes.object
}

export default StatisticsGrid
