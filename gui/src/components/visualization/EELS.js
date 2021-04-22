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
import React, {useState, useEffect, useMemo} from 'react'
import { useRecoilValue } from 'recoil'
import PropTypes from 'prop-types'
import { makeStyles, useTheme } from '@material-ui/core/styles'
import clsx from 'clsx'
import {
  Box
} from '@material-ui/core'
import Plot from '../visualization/Plot'
import { convertSI, convertSILabel, mergeObjects } from '../../utils'
import { withErrorHandler } from '../ErrorHandler'

const useStyles = makeStyles({
  root: {}
})

/**
 * Graph for EELS (electron energy loss specroscopy) data.
 */
function EELS({data, layout, aspectRatio, className, classes, unitsState, ...other}) {
  const [finalData, setFinalData] = useState(undefined)
  const units = useRecoilValue(unitsState)

  // Merge custom layout with default layout
  const tmpLayout = useMemo(() => {
    let defaultLayout = {
      yaxis: {
        title: {
          text: 'Electron count'
        }
      },
      xaxis: {
        showexponent: 'first',
        title: {
          text: `Electron energy loss (${convertSILabel('joule', units)})`
        }
      }
    }
    return mergeObjects(layout, defaultLayout)
  }, [layout, units])

  // Styles
  const style = useStyles(classes)
  const theme = useTheme()

  // The plotted data is loaded only after the first render as a side effect to
  // avoid freezing the UI
  useEffect(() => {
    if (data === undefined) {
      return
    }
    const plotData = []
    let energies = convertSI(data.energy, 'joule', units, false)
    plotData.push(
      {
        x: energies,
        y: data.count,
        type: 'scatter',
        mode: 'lines',
        line: {
          color: theme.palette.primary.main,
          width: 2
        }
      }
    )
    setFinalData(plotData)
  }, [data, theme.palette.primary.main, theme.palette.secondary.main, units])

  return (
    <Box className={clsx(style.root, className)}>
      <Plot
        data={finalData}
        layout={tmpLayout}
        aspectRatio={aspectRatio}
        floatTitle="EELS"
        fixedMargins={true}
        {...other}
      >
      </Plot>
    </Box>
  )
}

EELS.propTypes = {
  data: PropTypes.shape({
    count: PropTypes.arrayOf(PropTypes.number),
    energy: PropTypes.arrayOf(PropTypes.number)
  }),
  layout: PropTypes.object,
  aspectRatio: PropTypes.number,
  classes: PropTypes.object,
  className: PropTypes.string,
  unitsState: PropTypes.object // Recoil atom containing the unit configuration
}

export default withErrorHandler(EELS, 'Could not load EELS data.')
