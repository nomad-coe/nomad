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
import React, { createContext } from 'react'
import { makeStyles, useTheme } from '@material-ui/core/styles'
// import { Divider } from '@material-ui/core'
import PropTypes from 'prop-types'
import clsx from 'clsx'
import InputHeader from './InputHeader'
import InputTooltip from './InputTooltip'
import { InputGrid } from './InputGrid'
import { useSearchContext } from '../SearchContext'

/**
 * InputSection can be used to group together quantities that should be searched
 * together as an ES nested query:
 * https://www.elastic.co/guide/en/elasticsearch/reference/current/query-dsl-nested-query.html
 *
 * By wrapping search components in InputSection, the API calls will be
 * automatically performed as a nested query, and the visuals will change to
 * indicate the grouping.
 */
export const inputSectionContext = createContext()
const useStyles = makeStyles(theme => ({
  root: {
    width: '100%'
  },
  grid: {
    marginTop: theme.spacing(1)
  }
}))
const InputNestedObject = React.memo(({
  title,
  path,
  description,
  className,
  showHeader,
  classes,
  children,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useStyles({classes: classes, theme: theme})
  const { filterData } = useSearchContext()

  // Determine the description and units
  const def = filterData[path]
  const nested = def?.nested
  const repeats = def?.repeats
  const descFinal = description || def?.description || ''
  const labelFinal = title || def?.label

  return <InputTooltip>
    <div className={clsx(className, styles.root)} data-testid={testID}>
      {showHeader && <InputHeader
        quantity={path}
        label={labelFinal}
        description={descFinal}
        disableWidget
        disableStatistics
        />
      }
      <inputSectionContext.Provider value={{
        section: path,
        nested: nested,
        repeats: repeats
      }}>
        <InputGrid disablePadding className={styles.grid}>{children}</InputGrid>
      </inputSectionContext.Provider>
    </div>
  </InputTooltip>
})

InputNestedObject.propTypes = {
  title: PropTypes.string,
  path: PropTypes.string.isRequired,
  description: PropTypes.string,
  showHeader: PropTypes.bool,
  className: PropTypes.string,
  classes: PropTypes.object,
  children: PropTypes.node,
  'data-testid': PropTypes.string
}
InputNestedObject.defaultProps = {
  showHeader: true
}

export default InputNestedObject
