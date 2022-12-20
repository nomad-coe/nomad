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
import PropTypes from 'prop-types'
import clsx from 'clsx'
import InputHeader from './InputHeader'
import InputTooltip from './InputTooltip'
import { useSearchContext } from '../SearchContext'

export const inputSectionContext = createContext()

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    height: '100%',
    display: 'flex',
    justifyContent: 'center',
    flexDirection: 'column',
    boxSizing: 'border-box'
  },
  label: {
    marginBottom: `${theme.spacing(-0.5)}px !important`
  }
}))
const InputSection = React.memo(({
  label,
  section,
  description,
  className,
  disableHeader,
  classes,
  children,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useStyles({classes: classes, theme: theme})
  const { filterData } = useSearchContext()

  // Determine the description and units
  const def = filterData[section]
  const nested = def?.nested
  const repeats = def?.repeats
  const descFinal = description || def?.description || ''
  const labelFinal = label || def?.label

  return <InputTooltip>
    <div className={clsx(className, styles.root)} data-testid={testID}>
      {!disableHeader && <InputHeader
        quantity={section}
        label={labelFinal}
        description={descFinal}
        className={styles.label}
        disableWidget
        disableStatistics
      />}
      <inputSectionContext.Provider value={{
        section: section,
        nested: nested,
        repeats: repeats
      }}>
        {children}
      </inputSectionContext.Provider>
    </div>
  </InputTooltip>
})

InputSection.propTypes = {
  label: PropTypes.string,
  section: PropTypes.string.isRequired,
  description: PropTypes.string,
  disableHeader: PropTypes.bool,
  className: PropTypes.string,
  classes: PropTypes.object,
  children: PropTypes.node,
  'data-testid': PropTypes.string
}

export default InputSection
