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
import PropTypes from 'prop-types'
import { FilterSubMenu } from './FilterMenu'
import { InputGrid, InputGridItem } from '../input/InputGrid'
import InputRadio from '../input/InputRadio'
import { useApi } from '../../api'

const FilterSubMenuAccess = React.memo(({
  value,
  ...rest
}) => {
  const {api} = useApi()
  const authenticated = api?.keycloak?.authenticated

  return <FilterSubMenu value={value} {...rest}>
    <InputGrid>
      <InputGridItem xs={12}>
        <InputRadio
          quantity="visibility"
          label="Visibility"
          description="The visibility of the entry."
          initialValue={authenticated ? 'visible' : 'public'}
          options={{
            all: {label: 'All', disabled: !authenticated, tooltip: 'Consider all entries.'},
            public: {label: 'Public', disabled: false, tooltip: 'Consider all entries that can be publically downloaded, i.e. only published entries without embargo.'},
            visible: {label: 'Visible', disabled: !authenticated, tooltip: 'Consider all entries that are visible to you. This includes entries with embargo or unpublished entries that belong to you or are shared with you.'},
            shared: {label: 'Shared', disabled: !authenticated, tooltip: 'Only consider entries that belong to you or are shared with you.'},
            user: {label: 'User', disabled: !authenticated, tooltip: 'Only consider entries that belong to you.'},
            staging: {label: 'Unpublished', disabled: !authenticated, tooltip: 'Only search through unpublished entries.'}
          }}
        ></InputRadio>
      </InputGridItem>
    </InputGrid>
  </FilterSubMenu>
})
FilterSubMenuAccess.propTypes = {
  value: PropTypes.string
}

export default FilterSubMenuAccess
