
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
import React, { useEffect, useRef, useMemo } from 'react'
import { Paper } from '@material-ui/core'
import { makeStyles } from '@material-ui/core/styles'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import { useResizeDetector } from 'react-resize-detector'
import { isEmpty } from 'lodash'
import { MuuriComponent, getResponsiveStyle } from 'muuri-react'
import { useSearchContext } from '../SearchContext'

/**
 * Displays a summary of the anchored filter statistics in a grid.
 */
const widthMapping = {
  large: {
    'sm': 12,
    'md': 9,
    'lg': 8,
    'xl': 6
  },
  medium: {
    'sm': 6,
    'md': 6,
    'lg': 4,
    'xl': 3
  },
  small: {
    'sm': 6,
    'md': 3,
    'lg': 4,
    'xl': 3
  }
}
const useStyles = makeStyles(theme => {
  return {
    root: {},
    container: {
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
  const { useStatisticsValue, filterData } = useSearchContext()
  const { width, ref } = useResizeDetector()
  const statistics = useStatisticsValue()
  const gridRef = useRef()
  let size
  if (width > 1450) {
    size = 'xl'
  } else if (width > 1250) {
    size = 'lg'
  } else if (width > 950) {
    size = 'md'
  } else {
    size = 'sm'
  }

  // When the container size changes, the layout needs to be updated.
  useEffect(() => {
    if (gridRef.current) {
      gridRef.current.refreshItems()
      gridRef.current.layout()
    }
  }, [width])

  // Memoize the grid so that it gets rendered only when the contents or size
  // changes.
  const content = useMemo(() => {
    return (!isEmpty(statistics))
      ? <div className={styles.container}>
        <MuuriComponent
          ref={gridRef}
          dragEnabled
          layoutOnResize={false}
          dragHandle=".dragHandle"
          showDuration={0}
          hideDuration={0}
        >
          {Object.keys(statistics).map((filter) => {
            const config = filterData[filter].stats
            const layout = config.layout
            const muuriOuterItem = getResponsiveStyle({
              columns: widthMapping[layout.width][size] / 12,
              ratio: layout.ratio
            })
            return <div key={filter} style={muuriOuterItem}>
              <Paper className={styles.muuriInnerItem}>
                <config.component
                  quantity={filter}
                  visible
                  draggable
                  aggId="statistics"
                />
              </Paper>
            </div>
          })}
        </MuuriComponent>
      </div>
      : null
  }, [statistics, filterData, styles, size])

  return <div ref={ref} className={clsx(className, styles.root)}>
    {content}
  </div>
})

StatisticsGrid.propTypes = {
  className: PropTypes.string,
  classes: PropTypes.object
}

export default StatisticsGrid
