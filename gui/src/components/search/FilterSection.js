
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
import { makeStyles, useTheme } from '@material-ui/core/styles'
import { Typography } from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import RangeQuery from './RangeQuery'
import SelectQuery from './SelectQuery'
import searchQuantities from '../../searchQuantities'

/*
 * Used to compose a search component for quantities within an metainfo section.
 */
const useStaticStyles = makeStyles(theme => ({
  root: {
    padding: theme.spacing(2)
  },
  row: {
    display: 'flex',
    alignItems: 'center',
    margin: `${theme.spacing(1)}px 0`
  },
  title: {
  },
  contents: {
    padding: `${theme.spacing(0.5)}px ${theme.spacing(1)}px`
  }
}))
const SectionQuery = React.memo(({
  label,
  section,
  quantities,
  className,
  classes,
  units,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useStaticStyles({classes: classes, theme: theme})

  // Determine the displayed properties
  const content = Object.keys(quantities).map((key, index) => {
    const defCustom = quantities[key]
    const def = searchQuantities[`${section}.${key}`]
    const label = defCustom.label || key
    const options = defCustom.options || def?.type?.type_data
    if (defCustom.type === 'range') {
      return <RangeQuery key={key} label={label} quantity={`${section}.${key}`} units={units}/>
    } else if (defCustom.type === 'select') {
      return <SelectQuery key={key} label={label} quantity={`${section}.${key}`} options={options}/>
    }
  })

  return <div className={clsx(className, styles.root)} data-testid={testID}>
    <Typography variant="body1" className={styles.title}>
      {label}
    </Typography>
    <div className={styles.contents}>
      {content}
    </div>
  </div>
})

SectionQuery.propTypes = {
  /* Display name */
  label: PropTypes.string,
  /* Section path in metainfo */
  section: PropTypes.string,
  /* A mapping from section quantity names into the search component displayed
   * for them. */
  quantities: PropTypes.object,
  className: PropTypes.string,
  classes: PropTypes.object,
  units: PropTypes.object,
  'data-testid': PropTypes.string
}

export default SectionQuery
