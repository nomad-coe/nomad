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
import React, { useContext } from 'react'
import PropTypes from 'prop-types'
import { makeStyles } from '@material-ui/core/styles'
import { FilterSubMenu, filterMenuContext } from './FilterMenu'
import { InputGrid, InputGridItem } from '../input/InputGrid'
import InputPeriodicTable from '../input/InputPeriodicTable'
import InputField from '../input/InputField'
import InputRange from '../input/InputRange'

const useStyles = makeStyles(theme => ({
  grid: {
    marginTop: theme.spacing(2)
  },
  periodicTable: {
    height: '30rem'
  }
}))

const FilterSubMenuElements = React.memo(({
  id,
  ...rest
}) => {
  const {selected, open} = useContext(filterMenuContext)
  const visible = open && id === selected
  const styles = useStyles()

  return <FilterSubMenu id={id} {...rest}>
    <InputGrid className={styles.grid}>
      <InputGridItem xs={12}>
        <InputPeriodicTable
          quantity="results.material.elements"
          visible={visible}
          className={styles.periodicTable}
        />
      </InputGridItem>
      <InputGridItem xs={6}>
        <InputField
          quantity="results.material.chemical_formula_hill"
          visible={visible}
          disableOptions
        />
      </InputGridItem>
      <InputGridItem xs={6}>
        <InputField
          quantity="results.material.chemical_formula_iupac"
          visible={visible}
          disableOptions
        />
      </InputGridItem>
      <InputGridItem xs={6}>
        <InputField
          quantity="results.material.chemical_formula_reduced"
          visible={visible}
          disableOptions
        />
      </InputGridItem>
      <InputGridItem xs={6}>
        <InputField
          quantity="results.material.chemical_formula_anonymous"
          visible={visible}
          placeholder="E.g. A2B, A3B2"
          disableOptions
        />
      </InputGridItem>
      <InputGridItem xs={12}>
        <InputRange
          quantity="results.material.n_elements"
          visible={visible}
          step={1}
        />
      </InputGridItem>
    </InputGrid>
  </FilterSubMenu>
})
FilterSubMenuElements.propTypes = {
  id: PropTypes.string
}

export default FilterSubMenuElements
