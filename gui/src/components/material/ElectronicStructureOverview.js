import React, { useEffect, useState, useCallback } from 'react'
import PropTypes from 'prop-types'
import {
  Box,
  Card,
  CardHeader,
  CardContent
} from '@material-ui/core'
import DOS from '../visualization/DOS'
import BandStructure from '../visualization/BandStructure'
import BrillouinZone from '../visualization/BrillouinZone'
import { makeStyles } from '@material-ui/core/styles'
import { withApi } from '../api'

function ElectronicStructureOverview({data, range, className, classes, api, raiseError}) {
  const [dos, setDos] = useState()
  const [dosLayout, setDosLayout] = useState({yaxis: {range: range}})
  const [bs, setBs] = useState()
  const [bsLayout, setBsLayout] = useState({yaxis: {range: range}})

  // Styles
  const useStyles = makeStyles((theme) => {
    return {
      row: {
        display: 'flex',
        flexDirection: 'row',
        width: '100%'
      },
      bz: {
        flex: '1 1 25%'
      },
      bands: {
        flex: '1 1 50%'
      },
      dos: {
        flex: '1 1 25%'
      }
    }
  })
  const style = useStyles(classes)

  // Load the data parallelly from API on first render
  useEffect(() => {
    if (data === undefined) {
      return
    }
    // Check what data is available and request each in parallel
    let representatives = data.calculations.representatives
    let promises = []
    let requestedProperties = ['electronic_dos', 'electronic_band_structure']
    let availableProperties = []
    for (let property of requestedProperties) {
      if (representatives.hasOwnProperty(property)) {
        promises.push(
          api.encyclopediaCalculation(
            data.basic.material_id,
            representatives[property],
            {'properties': [property]}
          )
        )
        availableProperties.push(property)
      }
    }
    Promise.allSettled(promises).then((results) => {
      for (let i = 0; i < availableProperties.length; ++i) {
        let property = availableProperties[i]
        let result = results[i].value[property]
        if (property === 'electronic_dos') {
          setDos(result)
        }
        if (property === 'electronic_band_structure') {
          setBs(result)
        }
      }
    }).catch(error => {
      console.log(error)
      if (error.name === 'DoesNotExist') {
        raiseError(error)
      }
    })
  }, [data, api, raiseError])

  // Synchronize panning between BS/DOS plots
  const handleBSRelayouting = useCallback((event) => {
    let update = {
      yaxis: {
        range: [event['yaxis.range[0]'], event['yaxis.range[1]']]
      }
    }
    setDosLayout(update)
  }, [])
  const handleDOSRelayouting = useCallback((event) => {
    let update = {
      yaxis: {
        range: [event['yaxis.range[0]'], event['yaxis.range[1]']]
      }
    }
    setBsLayout(update)
  }, [])

  return (
    <Card>
      <CardHeader title="Electronic structure" />
      <CardContent>
        <Box className={style.row}>
          <BrillouinZone data={bs} className={style.bz} aspectRatio={1 / 2}></BrillouinZone>
          <BandStructure
            data={bs}
            layout={bsLayout}
            className={style.bands}
            aspectRatio={1}
            onRelayouting={handleBSRelayouting}
          ></BandStructure>
          <DOS
            data={dos}
            layout={dosLayout}
            className={style.dos}
            aspectRatio={1 / 2}
            onRelayouting={handleDOSRelayouting}
          ></DOS>
        </Box>
      </CardContent>
    </Card>
  )
}

ElectronicStructureOverview.propTypes = {
  data: PropTypes.object,
  range: PropTypes.array,
  className: PropTypes.string,
  classes: PropTypes.object,
  api: PropTypes.object,
  raiseError: PropTypes.func
}
ElectronicStructureOverview.defaultProps = {
  range: [-10, 20]
}

export default withApi(false, true)(ElectronicStructureOverview)
