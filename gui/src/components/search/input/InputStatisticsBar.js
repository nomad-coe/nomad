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
import { makeStyles } from '@material-ui/core/styles'
import PropTypes from 'prop-types'
const useStyles = makeStyles(theme => ({
  root: {
    width: '100%',
    height: '100%'
  },
  rectangle: {
    backgroundColor: theme.palette.secondary.main,
    width: '100%',
    height: '100%',
    '-webkit-transform': 'none',
    transform: 'none',
    transition: 'transform 250ms',
    willChange: 'transform'
  }
}))
const InputStatisticsBar = React.memo(({
  max,
  value,
  className,
  classes,
  'data-testid': testID
}) => {
  const styles = useStyles()

  return <div>
    <div className={styles.rectangle} styles={{transform: `scaleX(${value / max})`}}></div>
  </div>
})

InputStatisticsBar.propTypes = {
  max: PropTypes.number,
  value: PropTypes.number,
  className: PropTypes.string,
  classes: PropTypes.object,
  'data-testid': PropTypes.string
}

export default InputStatisticsBar
