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
import { Grid } from '@material-ui/core'
import { FilterSubMenu, filterMenuContext } from './FilterMenu'
import InputText from '../input/InputText'
import InputCheckboxes from '../input/InputCheckboxes'
import InputDateRange from '../input/InputDateRange'

const FilterSubMenuAuthor = React.memo(({
  value,
  ...rest
}) => {
  const {selected} = useContext(filterMenuContext)
  const visible = value === selected

  return <FilterSubMenu value={value} {...rest}>
    <Grid container spacing={2}>
      <Grid item xs={12}>
        <InputText
          label="author name"
          quantity="authors.name"
          visible={visible}
        />
      </Grid>
      <Grid item xs={12}>
        <InputCheckboxes
          label="external database"
          quantity="external_db"
          visible={visible}
        />
      </Grid>
      <Grid item xs={12}>
        <InputDateRange
          quantity="upload_create_time"
          visible={visible}
        />
      </Grid>
    </Grid>
  </FilterSubMenu>
})
FilterSubMenuAuthor.propTypes = {
  value: PropTypes.string
}

export default FilterSubMenuAuthor
