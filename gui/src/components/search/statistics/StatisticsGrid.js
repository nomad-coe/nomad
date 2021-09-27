
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
import useMediaQuery from '@material-ui/core/useMediaQuery'
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
      margin: -theme.spacing(1)
    },
    muuriInnerItem: {
      backgroundColor: theme.palette.background.paper,
      position: 'absolute',
      left: theme.spacing(1),
      top: theme.spacing(1),
      right: theme.spacing(1),
      bottom: theme.spacing(1),
      padding: theme.spacing(1),
      paddingLeft: theme.spacing(1.5),
      paddingRight: theme.spacing(1.5)
    }
  }
})
const StatisticsGrid = React.memo(({
  className,
  classes
}) => {
  const styles = useStyles(classes)
  const isXl = useMediaQuery('(min-width:1800px)')
  const isLg = useMediaQuery('(min-width:1600px)')
  const isMd = useMediaQuery('(min-width:1450px)')
  const size = isXl ? 'xl' : (isLg ? 'lg' : (isMd ? 'md' : 'sm'))
  const anchors = useAnchorsValue()

  return anchors.length > 0
    ? <div className={clsx(className, styles.root)}>
      <MuuriComponent
        dragEnabled
        dragHandle=".dragHandle"
        showDuration={0}
        hideDuration={0}
      >
        {anchors.map((filter) => {
          const config = filterData[filter].statConfig
          const layout = config.layout
          const muuriOuterItem = getResponsiveStyle({
            columns: widthMapping[layout.widthDefault][size] / 12,
            ratio: layout.ratioDefault
          })
          return <div key={filter} style={muuriOuterItem}>
            <Paper className={styles.muuriInnerItem}>
              <config.component
                quantity={filter}
                visible
                draggable
              />
            </Paper>
          </div>
        })}
      </MuuriComponent>
    </div>
    : null
})

StatisticsGrid.propTypes = {
  className: PropTypes.string,
  classes: PropTypes.object
}

export default StatisticsGrid
