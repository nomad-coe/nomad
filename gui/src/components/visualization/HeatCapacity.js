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
import React, {useEffect, useState, useMemo} from 'react'
import PropTypes from 'prop-types'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import { Box } from '@material-ui/core'
import clsx from 'clsx'
import Plot from '../visualization/Plot'
import { mergeObjects } from '../../utils'
import { convertSI, convertSILabel } from '../../units'
import { withErrorHandler } from '../ErrorHandler'

const useStyles = makeStyles({
  root: {}
})

const HeatCapacity = React.memo(({
  data,
  layout,
  aspectRatio,
  className,
  classes,
  placeholderStyle,
  units,
  'data-testid': testID,
  ...other
}) => {
  // Merge custom layout with default layout
  const finalLayout = useMemo(() => {
    let defaultLayout = {
      xaxis: {
        title: {
          text: `Temperature (${convertSILabel('kelvin', units)})`
        },
        zeroline: false
      },
      yaxis: {
        title: {
          text: `Heat capacity (${convertSILabel('joule/kelvin', units)})`
        },
        zeroline: false
      }
    }
    return mergeObjects(layout, defaultLayout)
  }, [layout, units])

  const [finalData, setFinalData] = useState(data === false ? data : undefined)
  const styles = useStyles(classes)
  const theme = useTheme()

  // Side effect that runs when the data that is displayed should change. By
  // running all this heavy stuff within useEffect (instead of e.g. useMemo),
  // the first render containing the placeholders etc. can be done as fast as
  // possible.
  useEffect(() => {
    if (!data) {
      return
    }

    // Convert units and determine range
    const temperatures = convertSI(data.temperatures, 'kelvin', units, false)
    const heatCapacities = convertSI(data.heat_capacities, 'joule/kelvin', units, false)

    // Create the final data that will be plotted.
    const plotData = [{
      x: temperatures,
      y: heatCapacities,
      type: 'scatter',
      mode: 'lines',
      line: {
        color: theme.palette.primary.main,
        width: 2
      }
    }]

    setFinalData(plotData)
  }, [data, units, theme])

  return (
    <Box className={clsx(styles.root, className)}>
      <Plot
        data={finalData}
        layout={finalLayout}
        aspectRatio={aspectRatio}
        floatTitle="Heat capacity"
        metaInfoLink={data?.m_path}
        data-testid={testID}
        {...other}
      >
      </Plot>
    </Box>
  )
})

HeatCapacity.propTypes = {
  data: PropTypes.oneOfType([
    PropTypes.bool, // Set to false to show NoData component
    PropTypes.shape({
      heat_capacities: PropTypes.array.isRequired,
      temperatures: PropTypes.array.isRequired,
      m_path: PropTypes.string // Path of the section containing the data in the Archive
    })
  ]),
  layout: PropTypes.object,
  aspectRatio: PropTypes.number,
  classes: PropTypes.object,
  className: PropTypes.string,
  placeholderStyle: PropTypes.string,
  noDataStyle: PropTypes.string,
  units: PropTypes.object, // Contains the unit configuration
  'data-testid': PropTypes.string
}
HeatCapacity.defaultProps = {
  type: 'electronic'
}

export default withErrorHandler(HeatCapacity, 'Could not load heat capacity.')
