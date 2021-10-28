
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
import searchQuantities from '../../../searchQuantities'
import InputHeader from './InputHeader'
import InputTooltip from './InputTooltip'

export const inputSectionContext = createContext()

const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    display: 'flex',
    justifyContent: 'center',
    flexDirection: 'column',
    boxSizing: 'border-box'
  },
  label: {
    marginBottom: theme.spacing(-0.5)
  }
}))
const InputSection = React.memo(({
  label,
  section,
  description,
  className,
  classes,
  children,
  'data-testid': testID
}) => {
  const theme = useTheme()
  const styles = useStyles({classes: classes, theme: theme})

  // Determine the description and units
  const def = searchQuantities[section]
  const desc = description || def?.description || ''
  const title = label || def?.name

  return <inputSectionContext.Provider value={{
    section: section
  }}>
    <InputTooltip>
      <div className={clsx(className, styles.root)} data-testid={testID}>
        <InputHeader
          quantity={section}
          label={title}
          description={desc}
          className={styles.label}
          variant="section"
          disableAggSize
          disableStatistics
        />
        {children}
      </div>
    </InputTooltip>
  </inputSectionContext.Provider>
})

InputSection.propTypes = {
  label: PropTypes.string,
  section: PropTypes.string.isRequired,
  description: PropTypes.string,
  className: PropTypes.string,
  classes: PropTypes.object,
  children: PropTypes.node,
  'data-testid': PropTypes.string
}

InputSection.defaultProps = {
  initialScale: 1
}

export default InputSection
