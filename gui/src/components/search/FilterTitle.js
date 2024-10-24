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
import React, { useMemo, useContext } from 'react'
import { Typography } from '@material-ui/core'
import PropTypes from 'prop-types'
import { useSearchContext } from './SearchContext'
import { inputSectionContext } from './input/InputNestedObject'
import { Unit } from '../units/Unit'
import { useUnitContext } from '../units/UnitContext'
import { DefinitionTitle } from '../DefinitionTitle'

/**
 * Title for a metainfo quantity or section that is used inside a search
 * context. By default the label, description and unit are automatically
 * retrieved from the filter config.
 */
const FilterTitle = React.memo(({
  quantity,
  label,
  description,
  unit,
  disableUnit,
  ...rest
}) => {
  const {filterData} = useSearchContext()
  const sectionContext = useContext(inputSectionContext)
  const {units} = useUnitContext()
  const section = sectionContext?.section

  // Create the final label
  const finalLabel = useMemo(() => {
    let finalLabel = label || filterData[quantity]?.label
    if (!disableUnit) {
      let finalUnit
      if (unit) {
        finalUnit = new Unit(unit).label()
      } else if (quantity && filterData[quantity]?.unit) {
        finalUnit = new Unit(filterData[quantity].unit).toSystem(units).label()
      }
      if (finalUnit) {
        finalLabel = `${finalLabel} (${finalUnit})`
      }
    }
    return finalLabel
  }, [filterData, quantity, units, label, unit, disableUnit])

  // Determine the final description
  let finalDescription = description || filterData[quantity]?.description || ''
  if (finalDescription && quantity) {
    finalDescription = (
      <>
        <Typography>{finalLabel}</Typography>
        <b>Description: </b>{finalDescription}<br/>
        <b>Path: </b>{quantity}
      </>
    )
  }

  return <DefinitionTitle
    label={finalLabel}
    description={finalDescription}
    section={section}
    {...rest}
  />
})

FilterTitle.propTypes = {
  quantity: PropTypes.string,
  label: PropTypes.string,
  description: PropTypes.string,
  unit: PropTypes.oneOfType([PropTypes.string, PropTypes.object]),
  disableUnit: PropTypes.bool
}

export default FilterTitle
