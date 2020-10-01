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

export default function DOS({data, layout, aspectRatio, className, classes, onRelayout, onAfterPlot, onRedraw, onRelayouting, unitsState}) {
  const [finalData, setFinalData] = useState(undefined)
  const units = useRecoilValue(unitsState)

  // Merge custom layout with default layout
  const tmpLayout = useMemo(() => {
    let defaultLayout = {
      yaxis: {
        title: {
          text: `Energy (${convertSILabel('joule', units)})`
        }
      }
    }
    return mergeObjects(layout, defaultLayout)
  }, [layout, units])

  // Styles
  const useStyles = makeStyles(
    {
      root: {
      }
    }
  )
  const style = useStyles(classes)
  const theme = useTheme()

  // The plotted data is loaded only after the first render as a side effect to
  // avoid freezing the UI
  useEffect(() => {
    if (data === undefined) {
      return
    }
    const norm = data.dos_energies_normalized === undefined ? '' : '_normalized'
    const energyName = 'dos_energies' + norm
    const valueName = 'dos_values' + norm
    const plotData = []
    if (data !== undefined) {
      let nChannels = data[valueName].length
      let energies = convertSI(data[energyName], 'joule', units, false)
      if (nChannels === 2) {
        plotData.push(
          {
            x: data[valueName][1],
            y: energies,
            type: 'scatter',
            mode: 'lines',
            line: {
              color: theme.palette.secondary.main,
              width: 2
            }
          }
        )
      }
      plotData.push(
        {
          x: data[valueName][0],
          y: energies,
          type: 'scatter',
          mode: 'lines',
          line: {
            color: theme.palette.primary.main,
            width: 2
          }
        }
      )
    }
    setFinalData(plotData)
  }, [data, theme.palette.primary.main, theme.palette.secondary.main, units])

  // Compute layout that depends on data.
  const computedLayout = useMemo(() => {
    if (data === undefined) {
      return {}
    }
    const norm = data.dos_energies_normalized !== undefined
    let defaultLayout = {
      xaxis: {
        title: {
          text: norm ? convertSILabel('states/joule/m^3/atom', units) : convertSILabel('states/joule/cell', units)
        }
      }
    }
    return defaultLayout
  }, [data, units])

  // Merge the given layout and layout computed from data
  const finalLayout = useMemo(() => {
    return mergeObjects(computedLayout, tmpLayout)
  }, [computedLayout, tmpLayout])

  return (
    <Box className={clsx(style.root, className)}>
      <Plot
        data={finalData}
        layout={finalLayout}
        aspectRatio={aspectRatio}
        floatTitle="Density of states"
        onRelayout={onRelayout}
        onAfterPlot={onAfterPlot}
        onRedraw={onRedraw}
        onRelayouting={onRelayouting}
      >
      </Plot>
    </Box>
  )
}

DOS.propTypes = {
  data: PropTypes.object, // section_dos
  layout: PropTypes.object,
  aspectRatio: PropTypes.number,
  classes: PropTypes.object,
  className: PropTypes.string,
  onAfterPlot: PropTypes.func,
  onRedraw: PropTypes.func,
  onRelayout: PropTypes.func,
  onRelayouting: PropTypes.func,
  unitsState: PropTypes.object // Recoil atom containing the unit configuration
}
