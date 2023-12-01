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
import React, { useState, useMemo, useEffect } from 'react'
import PropTypes from 'prop-types'
import assert from 'assert'
import { isNil, isArray } from 'lodash'
import { useTheme } from '@material-ui/core/styles'
import Plot from '../plotting/Plot'
import { QuantityTable, QuantityRow, QuantityCell } from '../Quantity'
import { ErrorHandler, withErrorHandler } from '../ErrorHandler'
import { diffTotal } from '../../utils'
import { Quantity } from '../units/Quantity'
import { Unit } from '../units/Unit'
import { useUnitContext } from '../units/UnitContext'
import { PropertyGrid, PropertyItem } from '../entry/properties/PropertyCard'

const energyUnit = new Unit('joule')
/**
 * Component for visualizing the results of a geometry optimization workflow.
 */
const GeometryOptimization = React.memo(({
  energies,
  convergence,
  className,
  classes,
  'data-testid': testID
}) => {
  const [finalData, setFinalData] = useState(!energies ? energies : undefined)
  const {units} = useUnitContext()
  const theme = useTheme()

  // Side effect that runs when the data that is displayed should change. By
  // running all this heavy stuff within useEffect (instead of e.g. useMemo),
  // the first render containing the placeholders etc. can be done as fast as
  // possible.
  useEffect(() => {
    if (!energies) {
      setFinalData(energies)
      return
    }

    // Convert energies into the correct units and calculate the total difference
    assert(
      isArray(energies),
      `Geometry optimization was expecting an array of energies, but instead got: ${energies}`
    )
    const energyDiffTotal = new Quantity(diffTotal(energies), energyUnit).toSystem(units).value()
    const convergenceCriteria = isNil(convergence?.convergence_tolerance_energy_difference)
      ? undefined
      : new Quantity(convergence.convergence_tolerance_energy_difference, energyUnit).toSystem(units).value()

    const steps = [...Array(energies.length).keys()]
    const energyDiff = []
    for (let i = 0; i < energyDiffTotal.length - 1; ++i) {
      energyDiff.push(Math.abs(energyDiffTotal[i + 1] - energyDiffTotal[i]))
    }
    const traces = [
      {
        x: steps,
        y: energyDiffTotal,
        name: 'Total change',
        type: 'scatter',
        showlegend: false,
        line: {
          color: theme.palette.primary.main,
          width: 2
        }
      },
      {
        x: steps.slice(1, steps.length),
        y: energyDiff,
        yaxis: 'y2',
        name: 'Abs. change per step',
        type: 'scatter',
        showlegend: false,
        line: {
          color: theme.palette.secondary.main,
          width: 2
        }
      }
    ]
    if (convergenceCriteria) {
      traces.push({
        x: [steps[0], steps[steps.length - 1]],
        y: [
          convergenceCriteria,
          convergenceCriteria
        ],
        yaxis: 'y2',
        name: 'Convergence criteria',
        showlegend: true,
        type: 'line',
        mode: 'lines',
        line: {
          color: theme.palette.secondary.main,
          width: 1,
          dash: '10px,10px'
        }
      })
    }
    setFinalData(traces)
  }, [energies, convergence, units, theme])

  const plotLayout = useMemo(() => {
    if (!energies) {
      return null
    }
    return {
      showlegend: true,
      legend: {
        x: 0,
        y: 1,
        xanchor: 'left'
      },
      xaxis: {
        showexponent: 'first',
        title: {
          text: 'Step number'
        },
        tickmode: 'auto',
        tickformat: ',d',
        autorange: false,
        range: [0, energies.length - 1],
        zeroline: false,
        showspikes: true,
        spikethickness: 2,
        spikedash: 'dot',
        spikecolor: '#999999',
        spikemode: 'across' },
      yaxis: {
        title: {
          text: `Total change (${energyUnit.toSystem(units).label()})`
        },
        tickfont: {
          color: theme.palette.primary.dark
        },
        autorange: true,
        zeroline: false
      },
      yaxis2: {
        title: {
          text: `Abs. change per step (${energyUnit.toSystem(units).label()})`
        },
        tickfont: {
          color: theme.palette.secondary.dark
        },
        type: 'log',
        autorange: true,
        zeroline: false,
        overlaying: 'y',
        side: 'right'
      }
    }
  }, [energies, theme, units])

  return <PropertyGrid>
    <PropertyItem title="Energy convergence" xs={12}>
      <ErrorHandler message='Could not load energies.'>
        <Plot
          data={finalData}
          layout={plotLayout}
          floatTitle="Energy convergence"
          data-testid={`${testID}-plot`}
        />
      </ErrorHandler>
    </PropertyItem>
    <PropertyItem title="Convergence results" xs={12} height="auto">
      <QuantityTable>
        <QuantityRow>
          <QuantityCell
            quantity='results.properties.geometry_optimization.final_energy_difference'
            value={convergence?.final_energy_difference}
          />
          <QuantityCell
            quantity='results.properties.geometry_optimization.final_displacement_maximum'
            value={convergence?.final_displacement_maximum}
          />
          <QuantityCell
            quantity='results.properties.geometry_optimization.final_force_maximum'
            value={convergence?.final_force_maximum}
          />
        </QuantityRow>
      </QuantityTable>
    </PropertyItem>
  </PropertyGrid>
})

GeometryOptimization.propTypes = {
  energies: PropTypes.oneOfType([
    PropTypes.bool,
    PropTypes.oneOf([undefined]),
    PropTypes.arrayOf(PropTypes.number)
  ]),
  convergence: PropTypes.oneOfType([
    PropTypes.bool,
    PropTypes.oneOf([undefined]),
    PropTypes.shape({
      final_energy_difference: PropTypes.number,
      final_displacement_maximum: PropTypes.number,
      final_force_maximum: PropTypes.number
    })
  ]),
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

GeometryOptimization.defaultProps = {
  'data-testid': 'geometry-optimization'
}

export default withErrorHandler('Could not load geometry optimization data.')(GeometryOptimization)
